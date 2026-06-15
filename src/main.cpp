#include <boost/program_options.hpp>
#include "Simulator.h"
#include "util_sim.h"
#include <map>

static std::string trim_main(std::string s)
{
    size_t first = s.find_first_not_of("\t\n\r ");
    if (first == std::string::npos)
        return "";
    size_t last = s.find_last_not_of("\t\n\r ");
    return s.substr(first, last - first + 1);
}

static std::vector<std::string> selected_init_lines(std::string init_states, std::string selection)
{
    std::set<std::string> selected;
    selection = trim_main(selection);
    bool select_all = selection == "" || selection == "all";
    if (!select_all)
    {
        std::stringstream selection_ss(selection);
        std::string label;
        while (getline(selection_ss, label, ','))
        {
            label = trim_main(label);
            if (label != "")
                selected.insert(label);
        }
    }

    std::vector<std::string> lines;
    std::set<std::string> found;
    std::stringstream init_ss(init_states);
    std::string line;
    while (getline(init_ss, line))
    {
        line = line.substr(0, line.find("//"));
        line = line.substr(0, line.find("#"));
        line = trim_main(line);
        if (line == "")
            continue;

        std::stringstream line_ss(line);
        std::string label;
        line_ss >> label;
        if (select_all || selected.find(label) != selected.end())
        {
            lines.push_back(line);
            found.insert(label);
        }
    }

    for (std::string label : selected)
        if (found.find(label) == found.end())
            std::cerr << "[warning]: Initial-state label '" << label << "' is not defined. It is skipped." << std::endl;

    return lines;
}

static std::string run_combined_init_simulation(
    int type, int shots, unsigned int seed, int r, bool isReorder, bool isQuery, bool isAlloc,
    const std::string& init_states, const std::string& selection, const std::string& obs_file,
    const std::string& qasm)
{
    std::ostringstream captured;
    std::streambuf *old_cout = std::cout.rdbuf(captured.rdbuf());
    {
        Simulator simulator(type, shots, seed, r, isReorder, isQuery, isAlloc);
        simulator.set_multi_initial_states(init_states);
        simulator.set_selected_initial_states(selection);
        if (obs_file != "")
            simulator.set_obs_file(obs_file);
        simulator.sim_qasm(qasm);
    }
    std::cout.rdbuf(old_cout);
    return captured.str();
}

static std::string run_sequential_init_simulation(
    int type, int shots, unsigned int seed, int r, bool isReorder, bool isQuery, bool isAlloc,
    const std::string& init_states, const std::string& selection, const std::string& obs_file,
    const std::string& qasm)
{
    std::ostringstream captured;
    std::streambuf *old_cout = std::cout.rdbuf(captured.rdbuf());
    {
        std::vector<std::string> lines = selected_init_lines(init_states, selection);
        for (std::string& line : lines)
        {
            Simulator seq_simulator(type, shots, seed, r, isReorder, isQuery, isAlloc);
            seq_simulator.set_multi_initial_states(line + "\n");
            if (obs_file != "")
                seq_simulator.set_obs_file(obs_file);
            seq_simulator.sim_qasm(qasm);
        }
    }
    std::cout.rdbuf(old_cout);
    return captured.str();
}

struct CountDistribution
{
    std::map<std::string, std::map<std::string, double>> probs;
};

static bool parse_quoted_value(const std::string& text, const std::string& key, size_t start, std::string *value)
{
    size_t key_pos = text.find(key, start);
    if (key_pos == std::string::npos)
        return false;
    size_t first_quote = text.find('"', key_pos + key.length());
    if (first_quote == std::string::npos)
        return false;
    size_t second_quote = text.find('"', first_quote + 1);
    if (second_quote == std::string::npos)
        return false;
    *value = text.substr(first_quote + 1, second_quote - first_quote - 1);
    return true;
}

static bool parse_counts_output(const std::string& output, CountDistribution *distribution)
{
    std::stringstream output_ss(output);
    std::string line;
    bool saw_counts = false;
    while (getline(output_ss, line))
    {
        if (line.find("\"counts\"") == std::string::npos)
            continue;

        std::string label;
        if (!parse_quoted_value(line, "\"initial_state\"", 0, &label))
            return false;

        size_t counts_pos = line.find("\"counts\"");
        size_t open = line.find('{', counts_pos);
        size_t close = line.find('}', open);
        if (open == std::string::npos || close == std::string::npos)
            return false;

        std::map<std::string, double> counts;
        double total = 0;
        size_t pos = open + 1;
        while (pos < close)
        {
            size_t key_start = line.find('"', pos);
            if (key_start == std::string::npos || key_start >= close)
                break;
            size_t key_end = line.find('"', key_start + 1);
            if (key_end == std::string::npos || key_end >= close)
                return false;
            std::string bitstring = line.substr(key_start + 1, key_end - key_start - 1);

            size_t colon = line.find(':', key_end);
            if (colon == std::string::npos || colon >= close)
                return false;
            size_t value_start = line.find_first_not_of(" \t", colon + 1);
            size_t value_end = line.find_first_of(",}", value_start);
            if (value_start == std::string::npos || value_end == std::string::npos)
                return false;
            double count = std::stod(line.substr(value_start, value_end - value_start));
            counts[bitstring] += count;
            total += count;
            pos = value_end + 1;
        }

        if (total <= 0)
            return false;
        for (auto& item : counts)
            distribution->probs[label][item.first] = item.second / total;
        saw_counts = true;
    }
    return saw_counts;
}

static bool compare_count_probabilities(
    const CountDistribution& combined,
    const CountDistribution& sequential,
    double tolerance,
    double *max_diff,
    std::string *max_diff_label,
    std::string *max_diff_output,
    bool *support_match,
    std::string *missing_from,
    std::string *missing_label,
    std::string *missing_output)
{
    *max_diff = 0;
    *support_match = true;
    std::set<std::string> labels;
    for (auto& item : combined.probs)
        labels.insert(item.first);
    for (auto& item : sequential.probs)
        labels.insert(item.first);

    bool pass = true;
    for (std::string label : labels)
    {
        std::set<std::string> outputs;
        auto combined_label = combined.probs.find(label);
        auto sequential_label = sequential.probs.find(label);
        if (combined_label == combined.probs.end())
        {
            if (*support_match)
            {
                *missing_from = "all-in-one";
                *missing_label = label;
                *missing_output = "";
            }
            *support_match = false;
            pass = false;
        }
        if (sequential_label == sequential.probs.end())
        {
            if (*support_match)
            {
                *missing_from = "sequential";
                *missing_label = label;
                *missing_output = "";
            }
            *support_match = false;
            pass = false;
        }
        if (combined_label != combined.probs.end())
            for (auto& item : combined_label->second)
                outputs.insert(item.first);
        if (sequential_label != sequential.probs.end())
            for (auto& item : sequential_label->second)
                outputs.insert(item.first);

        for (std::string output : outputs)
        {
            double combined_prob = 0;
            double sequential_prob = 0;
            bool in_combined = combined_label != combined.probs.end() && combined_label->second.find(output) != combined_label->second.end();
            bool in_sequential = sequential_label != sequential.probs.end() && sequential_label->second.find(output) != sequential_label->second.end();
            if (!in_combined)
            {
                if (*support_match)
                {
                    *missing_from = "all-in-one";
                    *missing_label = label;
                    *missing_output = output;
                }
                *support_match = false;
                pass = false;
            }
            if (!in_sequential)
            {
                if (*support_match)
                {
                    *missing_from = "sequential";
                    *missing_label = label;
                    *missing_output = output;
                }
                *support_match = false;
                pass = false;
            }

            if (in_combined)
                combined_prob = combined_label->second.at(output);
            if (in_sequential)
                sequential_prob = sequential_label->second.at(output);

            double diff = std::abs(combined_prob - sequential_prob);
            if (diff > *max_diff)
            {
                *max_diff = diff;
                *max_diff_label = label;
                *max_diff_output = output;
            }
            if (diff > tolerance)
                pass = false;
        }
    }
    return pass;
}

int main(int argc, char **argv)
{
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help",                                                     "produce help message")
    ("sim_qasm", po::value<std::string>()->implicit_value(""),   "simulate qasm file string")
    ("seed",     po::value<unsigned int>()->implicit_value(1),   "seed for random number generator")
    ("print_info",                                               "print simulation statistics such as runtime, memory, etc.")
    ("type",     po::value<unsigned int>()->default_value(0),    "the simulation type being executed.\n"
                                                                 "0: sampling mode (default option), where the sampled outcomes will be provided. \n"
                                                                 "1: all-amplitude mode, where the final state vector will be shown. \n"
                                                                 "2: query mode, where only the values of properties defined in obs_file will be provided.")
    ("shots",    po::value<unsigned int>()->default_value(1),    "the number of outcomes being sampled in \"sampling mode\". " )
    ("obs_file", po::value<std::string>(),                       "self-defined measurement operation file string (if any).")
    ("init_states", po::value<std::string>(),                     "multi-initial-state file. Each non-comment line is: <label> <bitstring> or <label> <amp0> <amp1> ...")
    ("select_init", po::value<std::string>()->default_value("all"), "initial-state label(s) to access after simulation. Use 'all' or comma-separated labels.")
    ("sequential_init",                                           "run selected initial states one by one instead of using selector variables. Useful for verification.")
    ("verify_init",                                               "run both all-in-one and sequential initial-state modes, then compare outputs.")
    ("verify_tol", po::value<double>()->default_value(1e-9),      "probability tolerance for --verify_init counts comparison.")
    ("dump_bdds", po::value<std::string>(),                       "write BDD Graphviz .dot dumps using the given filename prefix.")
    ("r",        po::value<unsigned int>()->default_value(32),   "integer bit size.")
    ("reorder",  po::value<bool>()->default_value(1),            "allow variable reordering or not.\n"
                                                                 "0: disable reordering.\n"
                                                                 "1: enable reordering (default option).")
    ("alloc",    po::value<bool>()->default_value(1),            "allocate new BDDs when overflow is detected.\n"
                                                                 "0: do not allocate new BDDs. This may lead to numerical errors.\n"
                                                                 "1: allocate new BDDs (default option).")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || !vm.count("sim_qasm") || argc == 1)
    {
	    std::cout << description << std::endl;
	    return 1;
	  }

    int type = vm["type"].as<unsigned int>(), shots = vm["shots"].as<unsigned int>(), r = vm["r"].as<unsigned int>();
    bool isReorder = vm["reorder"].as<bool>(), isAlloc = vm["alloc"].as<bool>();
    bool isQuery = bool(vm.count("obs_file"));

    std::random_device rd;
    unsigned int seed;
    if (vm.count("seed"))
        seed = vm["seed"].as<unsigned int>();
    else
        seed = rd();

    // start timer
    struct timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);

    assert(shots > 0);
    Simulator simulator(type, shots, seed, r, isReorder, isQuery, isAlloc);
    std::string init_states_str;
    std::string obs_file_str;

    if (vm.count("init_states"))
    {
        std::ifstream inFile;
        inFile.open(vm["init_states"].as<std::string>());
        if (!inFile)
        {
            std::cerr << "initial-state file doesn't exist." << std::endl;
            return -1;
        }
        std::stringstream initStream;
        initStream << inFile.rdbuf();
        init_states_str = initStream.str();
        simulator.set_multi_initial_states(init_states_str);
        simulator.set_selected_initial_states(vm["select_init"].as<std::string>());
    }

    if (vm.count("obs_file"))
    {
        std::ifstream inFile;
        inFile.open(vm["obs_file"].as<std::string>());
        if (!inFile)
        {
            std::cerr << "self-defined meaasurement operation file doesn't exist." << std::endl;
            return -1;
        }
        std::stringstream obsStream;
        obsStream << inFile.rdbuf();
        obs_file_str = obsStream.str();
        simulator.set_obs_file(obs_file_str);
    }

    if (vm.count("dump_bdds"))
        simulator.set_bdd_dump_prefix(vm["dump_bdds"].as<std::string>());

    // read in file into a string
    std::stringstream strStream;
    if (vm["sim_qasm"].as<std::string>() == "")
    {
        strStream << std::cin.rdbuf();    // read from std input; use Ctrl+D for ending
    }
    else
    {
        std::ifstream inFile;
        inFile.open(vm["sim_qasm"].as<std::string>()); //open the input file
        if (!inFile)
        {
            std::cerr << "qasm file doesn't exist." << std::endl;
            return -1;
        }
        strStream << inFile.rdbuf(); //read the file
    }
    std::string inFile_str = strStream.str(); //str holds the content of the file

    if (vm.count("verify_init"))
    {
        if (!vm.count("init_states"))
        {
            std::cerr << "--verify_init requires --init_states." << std::endl;
            return -1;
        }

        std::string selection = vm["select_init"].as<std::string>();
        std::string combined_output = run_combined_init_simulation(
            type, shots, seed, r, isReorder, isQuery, isAlloc,
            init_states_str, selection, obs_file_str, inFile_str);
        std::string sequential_output = run_sequential_init_simulation(
            type, shots, seed, r, isReorder, isQuery, isAlloc,
            init_states_str, selection, obs_file_str, inFile_str);

        CountDistribution combined_counts, sequential_counts;
        bool parsed_combined_counts = parse_counts_output(combined_output, &combined_counts);
        bool parsed_sequential_counts = parse_counts_output(sequential_output, &sequential_counts);
        bool compared_probabilities = parsed_combined_counts && parsed_sequential_counts;
        bool pass;
        double max_diff = 0;
        std::string max_diff_label;
        std::string max_diff_output;
        bool support_match = true;
        std::string missing_from;
        std::string missing_label;
        std::string missing_output;
        double verify_tol = vm["verify_tol"].as<double>();
        if (compared_probabilities)
        {
            pass = compare_count_probabilities(
                combined_counts, sequential_counts, verify_tol,
                &max_diff, &max_diff_label, &max_diff_output,
                &support_match, &missing_from, &missing_label, &missing_output);
        }
        else
        {
            pass = combined_output == sequential_output;
        }

        std::cout << "{ \"verification\": \"" << (pass ? "pass" : "fail") << "\"";
        if (compared_probabilities)
        {
            std::cout << ", \"comparison\": \"counts_probability\""
                      << ", \"tolerance\": " << verify_tol
                      << ", \"support_match\": " << (support_match ? "true" : "false")
                      << ", \"max_difference\": " << max_diff;
            if (!support_match)
            {
                std::cout << ", \"missing_from\": \"" << missing_from << "\""
                          << ", \"missing_label\": \"" << missing_label << "\"";
                if (missing_output != "")
                    std::cout << ", \"missing_output\": \"" << missing_output << "\"";
            }
            if (max_diff_label != "")
                std::cout << ", \"max_difference_label\": \"" << max_diff_label << "\""
                          << ", \"max_difference_output\": \"" << max_diff_output << "\"";
        }
        else
        {
            std::cout << ", \"comparison\": \"exact_output\"";
        }
        std::cout << " }" << std::endl;
        if (!pass)
        {
            std::cout << "[all-in-one output]" << std::endl;
            std::cout << combined_output;
            std::cout << "[sequential output]" << std::endl;
            std::cout << sequential_output;
        }
        return pass ? 0 : 2;
    }

    if (vm.count("sequential_init"))
    {
        if (!vm.count("init_states"))
        {
            std::cerr << "--sequential_init requires --init_states." << std::endl;
            return -1;
        }

        std::vector<std::string> lines = selected_init_lines(init_states_str, vm["select_init"].as<std::string>());
        for (std::string& line : lines)
        {
            Simulator seq_simulator(type, shots, seed, r, isReorder, isQuery, isAlloc);
            seq_simulator.set_multi_initial_states(line + "\n");
            if (obs_file_str != "")
                seq_simulator.set_obs_file(obs_file_str);
            seq_simulator.sim_qasm(inFile_str);
        }

        gettimeofday(&t2, NULL);
        elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
        elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

        double runtime = elapsedTime / 1000;
        size_t memPeak = getPeakRSS();
        if (vm.count("print_info"))
        {
            std::cout << "  Sequential verification mode" << std::endl;
            std::cout << "  #Initial states: " << lines.size() << std::endl;
            std::cout << "  Runtime: " << runtime << " seconds" << std::endl;
            std::cout << "  Peak memory usage: " << memPeak << " bytes" << std::endl;
        }
        return 0;
    }

    simulator.sim_qasm(inFile_str);

    //end timer
    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

    double runtime = elapsedTime / 1000;
    size_t memPeak = getPeakRSS();
    if (vm.count("print_info"))
        simulator.print_info(runtime, memPeak);

    return 0;
}
