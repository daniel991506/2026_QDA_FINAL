#include "Simulator.h"
#include "util_sim.h"

static std::string trim_copy(std::string s)
{
    size_t first = s.find_first_not_of("\t\n\r ");
    if (first == std::string::npos)
        return "";
    size_t last = s.find_last_not_of("\t\n\r ");
    return s.substr(first, last - first + 1);
}

void Simulator::set_multi_initial_states(std::string init_states)
{
    std::string line;
    std::stringstream inFile_ss(init_states);
    std::set<std::string> seen;

    while (getline(inFile_ss, line))
    {
        line = line.substr(0, line.find("//"));
        line = line.substr(0, line.find("#"));
        line = trim_copy(line);
        if (line == "")
            continue;

        std::stringstream line_ss(line);
        std::string label, entry;
        std::vector<std::string> entries;
        line_ss >> label;
        while (line_ss >> entry)
            entries.push_back(entry);
        if (label == "" || entries.empty())
        {
            std::cerr << "[warning]: Initial-state line '" << line << "' is ignored. Expected: <label> <bitstring> or <label> <amp0> <amp1> ..." << std::endl;
            continue;
        }
        if (label == "all")
        {
            std::cerr << "[warning]: Initial-state label 'all' is reserved. The line is ignored." << std::endl;
            continue;
        }
        if (seen.find(label) != seen.end())
        {
            std::cerr << "[warning]: Initial-state label '" << label << "' is duplicated. The later line is ignored." << std::endl;
            continue;
        }
        bool sparse_state = false;
        for (std::string& e : entries)
            if (e.find(':') != std::string::npos)
                sparse_state = true;
        bool basis_state = entries.size() == 1 && !sparse_state;
        bool legal = true;
        if (basis_state)
        {
            for (char c : entries[0])
            {
                if (c != '0' && c != '1')
                {
                    legal = false;
                    basis_state = false;
                    break;
                }
            }
        }
        if (!legal && entries.size() == 1 && !sparse_state)
        {
            std::cerr << "[warning]: Initial-state line '" << line << "' is ignored. Single-entry states must be bitstrings." << std::endl;
            continue;
        }
        seen.insert(label);
        multiInitLabels.push_back(label);
        multiInitEntries.push_back(entries);
        if (!basis_state || sparse_state)
            hasArbitraryInit = true;
    }

    if (multiInitLabels.empty())
        return;

    multiLabelBits = 0;
    int capacity = 1;
    while (capacity < multiInitLabels.size())
    {
        capacity *= 2;
        multiLabelBits++;
    }
    initScaleBits = std::max(0, std::min(24, r - 2));
}

void Simulator::set_selected_initial_states(std::string selection)
{
    selection = trim_copy(selection);
    if (selection == "" || selection == "all")
        return;

    std::stringstream selection_ss(selection);
    std::string label;
    while (getline(selection_ss, label, ','))
    {
        label = trim_copy(label);
        if (label != "")
            selectedInitLabels.push_back(label);
    }
}

void Simulator::set_obs_file(std::string obsfile)
{
    this->obsfile = obsfile;
}

void Simulator::set_bdd_dump_prefix(std::string prefix)
{
    bddDumpPrefix = prefix;
}

/**Function*************************************************************

  Synopsis    [Initailize simulator]

  Description [This function will set #qubits n, construct initial state, and enable dynamic reordering]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::init_simulator(int nQubits)
{
    n = nQubits; // set the number n here
    manager = Cudd_Init(n + multiLabelBits, n + multiLabelBits, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    measured_qubits_to_clbits = std::vector<std::vector<int>>(n, std::vector<int>(0));
    if (has_multi_initial_states())
    {
        for (int i = 0; i < multiInitEntries.size(); i++)
        {
            std::vector<std::string>& entries = multiInitEntries[i];
            bool sparse_state = false;
            for (std::string& e : entries)
                if (e.find(':') != std::string::npos)
                    sparse_state = true;
            if (entries.size() == 1 && !sparse_state && entries[0].length() != n)
            {
                std::cerr << "[Error]: Initial state '" << entries[0] << "' has length " << entries[0].length() << ", but qreg has " << n << " qubits." << std::endl;
                assert(entries[0].length() == n);
            }
            if (!sparse_state && entries.size() != 1)
            {
                if (n >= 63)
                {
                    std::cerr << "[Error]: Dense amplitude initial state '" << multiInitLabels[i] << "' is too large for " << n << " qubits. Use sparse entries like <bitstring>:<amplitude>." << std::endl;
                    assert(n < 63);
                }
                unsigned long long nEntries = 1ULL << n;
                if (entries.size() != nEntries)
                {
                    std::cerr << "[Error]: Initial state '" << multiInitLabels[i] << "' has " << entries.size() << " amplitudes, but qreg with " << n << " qubits requires " << nEntries << "." << std::endl;
                    assert(entries.size() == nEntries);
                }
            }
        }
        if (hasArbitraryInit)
        {
            k = 2 * initScaleBits;
        }
        init_multi_state();
        dump_all_bdds("initial_multi_input");
    }
    else
    {
        int *constants = new int[n];
        for (int i = 0; i < n; i++)
            constants[i] = 0; // TODO: costom initial state
        init_state(constants);
        delete[] constants;
        dump_all_bdds("initial_input");
    }
    if (isReorder) Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);
}


/**Function*************************************************************

  Synopsis    [parse and simulate the qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm_file(std::string qasm)
{
    std::string inStr;
    std::stringstream inFile_ss(qasm);
    while (getline(inFile_ss, inStr))
    {
        inStr = inStr.substr(0, inStr.find("//"));
        if (inStr.find_first_not_of("\t\n\r ") != std::string::npos)
        {        
            std::stringstream inStr_ss(inStr);
            getline(inStr_ss, inStr, ' ');
            if (inStr == "qreg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                init_simulator(stoi(inStr));
                
                if (sim_type == 1 && n > 50) 
                {
                    std::cerr << "[Error]: The all-amplitude mode will print the whole state vector, which is too long for large qubit number. Please consider using the sampling mode." << std::endl;
                    assert(sim_type != 1 || n <= 50); 
                }
            }
            else if (inStr == "creg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                nClbits = stoi(inStr);                
            }
            else if (inStr == "OPENQASM"){;}
            else if (inStr == "include"){;}
            else if (inStr == "measure")
            {
                isMeasure = 1;
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                int qIndex = stoi(inStr);
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                int cIndex = stoi(inStr);
                measure(qIndex, cIndex);
            }
            else
            {
                if (inStr == "x")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliX(stoi(inStr));
                }
                else if (inStr == "y")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliY(stoi(inStr));
                }
                else if (inStr == "z")
                {
                    std::vector<int> iqubit(1);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    iqubit[0] = stoi(inStr);
                    PauliZ(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "h")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    Hadamard(stoi(inStr));
                }
                else if (inStr == "s")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(2, iqubit);
                }
                else if (inStr == "sdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift_dagger(-2, iqubit);
                }
                else if (inStr == "t")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(4, iqubit);
                }
                else if (inStr == "tdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift_dagger(-4, iqubit);
                }
                else if (inStr == "rx(pi/2)")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    rx_pi_2(stoi(inStr));
                }
                else if (inStr == "ry(pi/2)")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    ry_pi_2(stoi(inStr));
                }
                else if (inStr == "cx")
                {
                    std::vector<int> cont(1);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    cont[0] = stoi(inStr);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    targ = stoi(inStr);
                    Toffoli(targ, cont, ncont);
                    cont.clear();
                    ncont.clear();
                }
                else if (inStr == "cz")
                {
                    std::vector<int> iqubit(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                    }
                    PauliZ(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "swap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(0);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont);
                    cont.clear();
                }
                else if (inStr == "cswap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(1);
                    for (int i = 0; i < 3; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            cont[i] = stoi(inStr);
                        else if (i == 1)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont);
                    cont.clear();
                }
                else if (inStr == "ccx" || inStr == "mcx")
                {
                    std::vector<int> cont(0);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    while(getline(inStr_ss, inStr, ']'))
                    {
                        cont.push_back(stoi(inStr));
                        getline(inStr_ss, inStr, '[');
                    }
                    targ = cont.back();
                    cont.pop_back();
                    Toffoli(targ, cont, ncont);
                    cont.clear();
                    ncont.clear();
                }
                else
                {
                    std::cerr << std::endl
                            // << "[warning]: Gate \'" << inStr << "\' is not supported in this simulator. The gate is ignored ..." << std::endl;
                            << "[warning]: Syntax \'" << inStr << "\' is not supported in this simulator. The line is ignored ..." << std::endl;
                }
            }
        }
    }
    if (isReorder) Cudd_AutodynDisable(manager);
}

/**Function*************************************************************

  Synopsis    [simulate the circuit described by a qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm(std::string qasm)
{
    sim_qasm_file(qasm); // simulate

    if (has_multi_initial_states())
    {
        if (sim_type == 0 && isMeasure == 0)
        {
            std::cerr << "[Error]: no measurement detected. Cannot do sampling." << std::endl;
            assert(sim_type != 0 || isMeasure != 0);
        }
        if (sim_type == 1 && isMeasure == 1 && shots != 1)
        {
            shots = 1;
            std::cerr << "[Warning]: shot number is limited to 1 in all_amplitude mode." << std::endl;
        }

        if (bddDumpPrefix != "")
        {
            dump_all_bdds("final_multi_state");
            create_bigBDD();
            dump_big_bdd("final_merged_all_in_one");
            Cudd_RecursiveDeref(manager, bigBDD);
            bigBDD = nullptr;
        }

        for (std::string label : labels_to_access())
        {
            if (!select_initial_state(label))
                continue;

            if (sim_type == 0)
            {
                create_bigBDD();
                dump_big_bdd("final_merged_" + label);
                measurement();
                print_results(label);
                if (obsfile != "")
                {
                    std::cout << "\"initial_state\": \"" << label << "\"" << std::endl;
                    measurement_obs(obsfile);
                }
            }
            else if (sim_type == 1)
            {
                if (isMeasure || isQuery)
                {
                    create_bigBDD();
                    dump_big_bdd("final_merged_" + label);
                    if (isMeasure)
                        measurement();
                }
                getStatevector();
                print_results(label);
                if (obsfile != "")
                {
                    std::cout << "\"initial_state\": \"" << label << "\"" << std::endl;
                    measurement_obs(obsfile);
                }
            }
            else if (sim_type == 2)
            {
                create_bigBDD();
                dump_big_bdd("final_merged_" + label);
                if (obsfile != "")
                {
                    std::cout << "\"initial_state\": \"" << label << "\"" << std::endl;
                    measurement_obs(obsfile);
                }
            }
            else
            {
                std::cerr << "[Error]: unknown sim_type." << std::endl;
                exit(-1);
            }

            reset_access_state();
            restore_multi_state();
        }
        return;
    }

    if (sim_type == 0 && isMeasure == 0)
    {
        std::cerr << "[Error]: no measurement detected. Cannot do sampling." << std::endl;
        assert(sim_type != 0 || isMeasure != 0);
    }
    if (sim_type == 1)
    {
        
        if (isMeasure == 1)
        {
            std::cerr << "[Warning]: measurement detected. The final statevector will collapse based on the measurement outcome." << std::endl;
            if (shots != 1)
            {
                shots = 1;
                std::cerr << "[Warning]: shot number is limited to 1 in all_amplitude mode." << std::endl;
            }
        }    
        else 
        {
            if (shots != 1)
            {
                std::cerr << "[Warning]: no measurement detected. The --shots argument is ignored." << std::endl;
            }
        }  
    }

    // measure based on simulator type
    if (sim_type == 0) // sampling mode
    {
        create_bigBDD();
        dump_big_bdd("final_merged");
        measurement();
        print_results();
        if (obsfile != "")
            measurement_obs(obsfile);
    }
    else if (sim_type == 1) // all_amplitude mode
    {
        if (isMeasure)
        {
            create_bigBDD();
            dump_big_bdd("final_merged");
            measurement();
        }
        else if (isQuery)
        {
            create_bigBDD();
            dump_big_bdd("final_merged");
        }
        else if (bddDumpPrefix != "")
        {
            create_bigBDD();
            dump_big_bdd("final_merged");
            Cudd_RecursiveDeref(manager, bigBDD);
            bigBDD = nullptr;
        }
        getStatevector();
        print_results();
        if (obsfile != "")
            measurement_obs(obsfile);
    }
    else if (sim_type == 2)
    {
        create_bigBDD();
        dump_big_bdd("final_merged");
        if (obsfile != "")
            measurement_obs(obsfile);
    }
    else
    {
        std::cerr << "[Error]: unknown sim_type." << std::endl;
        exit(-1);
    }
}



/**Function*************************************************************

  Synopsis    [print state vector and distribution of sampled outcomes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::print_results(std::string label)
{
    // write output string based on state_count and statevector
    std::unordered_map<std::string, int>::iterator it;
    
    run_output = "{";
    if (state_count.begin() != state_count.end()){
        run_output += "\"counts\": { ";
        for (it = state_count.begin(); it != state_count.end(); it++)
        {
            if (std::next(it) == state_count.end())
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second);
            else
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second) + ", ";
        }
        run_output += " }";
        run_output += (statevector != "null") ? ", " : ""; 
    }    

    run_output += (statevector != "null") ? "\"statevector\": " + statevector + " }" : " }";
    //return;
    if (label == "")
        std::cout << run_output << std::endl;
    else
        std::cout << "{\"initial_state\": \"" << label << "\", \"result\": " << run_output << "}" << std::endl;
}
