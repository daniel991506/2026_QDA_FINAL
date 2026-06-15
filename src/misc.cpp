#include "Simulator.h"
#include "util_sim.h"
#include <algorithm>
#include <cctype>

static bool is_bitstring(const std::string& value)
{
    if (value.empty())
        return false;
    for (char c : value)
        if (c != '0' && c != '1')
            return false;
    return true;
}

static std::complex<double> parse_complex_value(std::string raw, bool *ok)
{
    *ok = true;
    raw.erase(remove_if(raw.begin(), raw.end(), ::isspace), raw.end());
    if (raw.empty())
    {
        *ok = false;
        return {0, 0};
    }
    if (raw.front() == '(' && raw.back() == ')')
        raw = raw.substr(1, raw.length() - 2);

    size_t i_pos = raw.find('i');
    if (i_pos == std::string::npos)
    {
        char *end = nullptr;
        double real = strtod(raw.c_str(), &end);
        if (*end != '\0')
            *ok = false;
        return {real, 0};
    }
    if (i_pos != raw.length() - 1 || raw.find('i', i_pos + 1) != std::string::npos)
    {
        *ok = false;
        return {0, 0};
    }

    std::string without_i = raw.substr(0, raw.length() - 1);
    size_t split = std::string::npos;
    for (size_t i = 1; i < without_i.length(); i++)
    {
        if (without_i[i] == '+' || without_i[i] == '-')
            split = i;
    }

    std::string real_part = "0";
    std::string imag_part = without_i;
    if (split != std::string::npos)
    {
        real_part = without_i.substr(0, split);
        imag_part = without_i.substr(split);
    }
    if (imag_part == "" || imag_part == "+")
        imag_part = "1";
    else if (imag_part == "-")
        imag_part = "-1";

    char *end_real = nullptr;
    char *end_imag = nullptr;
    double real = strtod(real_part.c_str(), &end_real);
    double imag = strtod(imag_part.c_str(), &end_imag);
    if (*end_real != '\0' || *end_imag != '\0')
        *ok = false;
    return {real, imag};
}

static long long encode_fixed_point(double value, int scale_bits, int r, bool *ok)
{
    double scaled = value * pow(2.0, scale_bits);
    long long encoded = llround(scaled);
    long long min_value = -(1LL << (r - 1));
    long long max_value = (1LL << (r - 1)) - 1;
    if (encoded < min_value || encoded > max_value)
        *ok = false;
    return encoded;
}

static bool integer_bit(long long value, int bit)
{
    unsigned long long encoded = static_cast<unsigned long long>(value);
    return ((encoded >> bit) & 1ULL) != 0;
}

static DdNode* build_cube(DdManager *manager, int first_var, int bit_count, unsigned long long value)
{
    DdNode *cube = Cudd_ReadOne(manager);
    Cudd_Ref(cube);
    for (int j = bit_count - 1; j >= 0; j--)
    {
        bool bit = (value >> j) & 1ULL;
        DdNode *tmp = Cudd_bddAnd(manager, bit ? Cudd_bddIthVar(manager, first_var + j) : Cudd_Not(Cudd_bddIthVar(manager, first_var + j)), cube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, cube);
        cube = tmp;
    }
    return cube;
}

static DdNode* build_bitstring_cube(DdManager *manager, int first_var, const std::string& bits)
{
    DdNode *cube = Cudd_ReadOne(manager);
    Cudd_Ref(cube);
    for (int j = bits.length() - 1; j >= 0; j--)
    {
        bool bit = bits[bits.length() - 1 - j] == '1';
        DdNode *tmp = Cudd_bddAnd(manager, bit ? Cudd_bddIthVar(manager, first_var + j) : Cudd_Not(Cudd_bddIthVar(manager, first_var + j)), cube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, cube);
        cube = tmp;
    }
    return cube;
}

static bool has_sparse_amplitudes(const std::vector<std::string>& entries)
{
    for (const std::string& entry : entries)
        if (entry.find(':') != std::string::npos)
            return true;
    return false;
}


/**Function*************************************************************

  Synopsis    [initialize state vector by a basis state]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::init_state(int *constants)
{
    DdNode *var, *tmp;
    All_Bdd = new DdNode **[w];
    for (int i = 0; i < w; i++)
        All_Bdd[i] = new DdNode *[r];

    for (int i = 0; i < r; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < w - 1; j++)
            {
                All_Bdd[j][i] = Cudd_Not(Cudd_ReadOne(manager));
                Cudd_Ref(All_Bdd[j][i]);
            }
            All_Bdd[w - 1][i] = Cudd_ReadOne(manager);
            Cudd_Ref(All_Bdd[w - 1][i]);
            for (int j = n - 1; j >= 0; j--)
            {
                var = Cudd_bddIthVar(manager, j);
                if (constants[j] == 0)
                    tmp = Cudd_bddAnd(manager, Cudd_Not(var), All_Bdd[w - 1][i]);
                else
                    tmp = Cudd_bddAnd(manager, var, All_Bdd[w - 1][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, All_Bdd[w - 1][i]);
                All_Bdd[w - 1][i] = tmp;
            }
        }
        else
        {
            for (int j = 0; j < w; j++)
            {
                All_Bdd[j][i] = Cudd_Not(Cudd_ReadOne(manager));
                Cudd_Ref(All_Bdd[j][i]);
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [allocate new BDDs for each integer vector]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::alloc_BDD(DdNode ***Bdd, bool extend)
{
    DdNode *tmp;

    DdNode ***W = new DdNode **[w];
    for (int i = 0; i < w; i++)
        W[i] = new DdNode *[r];

    for (int i = 0; i < r - inc; i++)
        for (int j = 0; j < w; j++)
            W[j][i] = Bdd[j][i];

    for (int i = 0; i < w; i++)
        delete[] Bdd[i];

    for (int i = 0; i < w; i++)
        Bdd[i] = W[i];

    if (extend)
    {
        for (int i = r - inc; i < r; i++)
        {
            for (int j = 0; j < w; j++)
            {
                Bdd[j][i] = Cudd_ReadOne(manager);
                Cudd_Ref(Bdd[j][i]);
                tmp = Cudd_bddAnd(manager, Bdd[j][r - inc - 1], Bdd[j][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, Bdd[j][i]);
                Bdd[j][i] = tmp;
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [Drop LSB and shift right by 1 bit]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::dropLSB(DdNode ***Bdd)
{
    DdNode *tmp;

    for (int i = 0; i < w; i++)
    {
        Cudd_RecursiveDeref(manager, Bdd[i][0]); // drop LSB
        // right shift
        for (int j = 1; j < r; j++)
        {
            Bdd[i][j - 1] = Bdd[i][j];
        }
        // sign extension
        Bdd[i][r - 1] = Cudd_ReadOne(manager);
        Cudd_Ref(Bdd[i][r - 1]);
        tmp = Cudd_bddAnd(manager, Bdd[i][r - 2], Bdd[i][r - 1]);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, Bdd[i][r - 1]);
        Bdd[i][r - 1] = tmp;
    }
}

/**Function*************************************************************

  Synopsis    [detect overflow in integer vectors]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int Simulator::overflow3(DdNode *g, DdNode *h, DdNode *crin)
{
    DdNode *tmp, *dd1, *dd2;
    int overflow;

    dd1 = Cudd_bddXor(manager, g, crin);
    Cudd_Ref(dd1);

    dd2 = Cudd_bddXnor(manager, g, h);
    Cudd_Ref(dd2);

    tmp = Cudd_bddAnd(manager, dd1, dd2);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, dd1);
    Cudd_RecursiveDeref(manager, dd2);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(manager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [detect overflow in integer vectors -- for the case that h is 0]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int Simulator::overflow2(DdNode *g, DdNode *crin){
    DdNode *tmp;
    int overflow;

    tmp = Cudd_bddAnd(manager, Cudd_Not(g), crin);
    Cudd_Ref(tmp);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(manager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [decode and print each entry of the state vector]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::decode_entries()
{
    double oneroot2 = 1 / sqrt(2);
    double H_factor = pow(oneroot2, k);
    double re = 0, im = 0;
    int *assign = new int[n];
    unsigned long long nEntries = pow(2, n);
    int oneEntry;
    long long int_value = 0;
    DdNode *tmp;

    for (int i = 0; i < n; i++) //initialize assignment
        assign[i] = 0;

    std::cout << "Amplitudes of the Computational Basis States:" << std::endl;

    for (unsigned long long i = 0; i < nEntries; i++) // compute every entry
    {
        re = 0;
        im = 0;
        for (int j = 0; j < w; j++) // compute every complex value
        {
            int_value = 0;
            for (int h = 0; h < r; h++) // compute every integer
            {
                tmp = Cudd_Eval(manager, All_Bdd[j][h], assign);
                Cudd_Ref(tmp);
                oneEntry = !(Cudd_IsComplement(tmp));
                Cudd_RecursiveDeref(manager, tmp);
                if (h == r - 1)
                    int_value -= oneEntry * pow(2, h);
                else
                    int_value += oneEntry * pow(2, h);
            }
            /* translate to re and im */
            re += int_value * cos((double) (w - j - 1)/w * PI);
            im += int_value * sin((double) (w - j - 1)/w * PI);
        }
        re *= H_factor;
        im *= H_factor;
        std::cout << i << ": " << re << " + " << im << "i" << std::endl;
        full_adder_plus_1(n, assign);
    }

    delete[] assign;
}

/**Function*************************************************************

  Synopsis    [reorder BDDs]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::reorder()
{
    int reorder_signal = Cudd_ReduceHeap(manager, CUDD_REORDER_SYMM_SIFT, 0);
    if (!reorder_signal)
        std::cout << "reorder fails" << std::endl;
}

/**Function*************************************************************

  Synopsis    [update max #nodes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::nodecount()
{
    unsigned long NodeCount_new = Cudd_ReadNodeCount(manager);
    if (NodeCount_new > NodeCount)
         NodeCount = NodeCount_new;
}

/**Function*************************************************************

  Synopsis    [print statistics]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::print_info(double runtime, size_t memPeak)
{
    std::cout << "  Runtime: " << runtime << " seconds" << std::endl;
    std::cout << "  Peak memory usage: " << memPeak << " bytes" << std::endl; //unit in bytes
    std::cout << "  #Applied gates: " << gatecount << std::endl;
    std::cout << "  Max #nodes: " << NodeCount << std::endl;
    std::cout << "  Integer bit size: " << r << std::endl;
    std::cout << "  Accuracy loss: " << error << std::endl;
    // std::cout << "  #Integers: " << w << std::endl;

    // std::unordered_map<std::string, int>::iterator it;
    // std::cout << "  Measurement: " << std::endl;
    // for(it = state_count.begin(); it != state_count.end(); it++)
    //     std::cout << "      " << it->first << ": " << it->second << std::endl;
}

/**Function*************************************************************

  Synopsis    [parse an infix Boolean function to postfix]

  Description [use shunting yard algorithm]

  SideEffects []

  SeeAlso     []

***********************************************************************/

std::vector<std::string> Simulator::boolean_parser(std::string& inStr)
{
    inStr = std::regex_replace(inStr, std::regex("!"), " not ");
    inStr = std::regex_replace(inStr, std::regex("~"), " not ");
    inStr = std::regex_replace(inStr, std::regex("\\*"), " and ");
    inStr = std::regex_replace(inStr, std::regex("&"), " and ");
    inStr = std::regex_replace(inStr, std::regex("\\+"), " or ");
    inStr = std::regex_replace(inStr, std::regex("\\|"), " or ");
    inStr = std::regex_replace(inStr, std::regex("\\^"), " xor ");
    inStr = std::regex_replace(inStr, std::regex("\\("), " ( ");
    inStr = std::regex_replace(inStr, std::regex("\\)"), " ) ");
    
    std::vector<std::string> result;
    std::map<std::string, int> operator_priority = { {"not", 4},
                                                     {"xor", 3},
                                                     {"and", 2},
                                                     {"or",  1}, };
        
    std::stringstream inStr_ss(inStr);
    std::stack<std::string> waiting_operators;
    int nLeftTuple = 0;
    while(getline(inStr_ss, inStr, ' '))
    {
        if (inStr == "") continue;
                
        if (inStr == "(")
        {
            waiting_operators.push(inStr);
            nLeftTuple++;
        }
        else if (inStr == ")")
        {
            assert(nLeftTuple > 0);
            while (waiting_operators.top() != "(") 
            {
                result.emplace_back(waiting_operators.top());
                waiting_operators.pop();
            }
            waiting_operators.pop();
            nLeftTuple--;
        }
        else if (operator_priority.count(inStr) != 0)
        {
            while (!waiting_operators.empty() && waiting_operators.top() != "(" && operator_priority[inStr] <= operator_priority[waiting_operators.top()])
            {
                result.emplace_back(waiting_operators.top());
                waiting_operators.pop();
            }
            waiting_operators.push(inStr);
        }
        else if (inStr.find( "[" ) != std::string::npos) 
        {
            std::stringstream inStr_sss(inStr);
            getline(inStr_sss, inStr, '[');
            getline(inStr_sss, inStr, ']');
            assert(inStr != "");
            result.emplace_back(inStr);
        }
        else  // self-defined variables
        {
            result.emplace_back(inStr);
        }
    }
    assert(nLeftTuple == 0);
    
    while (!waiting_operators.empty())
    {
        result.emplace_back(waiting_operators.top());
        waiting_operators.pop();
    }
    return result;
}

/**Function*************************************************************

  Synopsis    [compare two node lists and return a node representing 
               whether two node lists represents the same integers ]

  Description []
               
  SideEffects [returned node already referenced]

  SeeAlso     []

***********************************************************************/
DdNode* Simulator::node_equiv(std::vector<DdNode*>& int_1, std::vector<DdNode*>& int_2)
{
    DdNode* result = Cudd_ReadOne(manager);
    Cudd_Ref(result);
    
    for (int i = 0; i < int_1.size(); ++i)
    {
        DdNode* tmp = Cudd_bddXnor(manager, int_1[i], int_2[i]);
        Cudd_Ref(tmp);
        DdNode* new_result = Cudd_bddAnd(manager, tmp, result);
        Cudd_Ref(new_result);
        Cudd_RecursiveDeref(manager, tmp);
        Cudd_RecursiveDeref(manager, result);
        result = new_result;
    }
    
    return result;
}

/**Function*************************************************************

  Synopsis    [compare two node lists and return a node representing 
               whether the first node list represents a larger integer ]

  Description []
               
  SideEffects [returned node already referenced]

  SeeAlso     []

***********************************************************************/
DdNode* Simulator::node_larger(std::vector<DdNode*>& int_1, std::vector<DdNode*>& int_2)
{
    DdNode* result = Cudd_Not(Cudd_ReadOne(manager));
    Cudd_Ref(result);
    
    DdNode* even = Cudd_ReadOne(manager);
    Cudd_Ref(even);
    
    for (int i = 0; i < int_1.size(); ++i)
    {
        DdNode* tmp1 = Cudd_bddAnd(manager, int_1[i], Cudd_Not(int_2[i]));
        Cudd_Ref(tmp1);
        DdNode* new_larger = Cudd_bddAnd(manager, tmp1, even);
        Cudd_Ref(new_larger);
        DdNode* new_result = Cudd_bddOr(manager, result, new_larger);
        Cudd_Ref(new_result);
        Cudd_RecursiveDeref(manager, tmp1);
        Cudd_RecursiveDeref(manager, new_larger);
        Cudd_RecursiveDeref(manager, result);
        result = new_result;
        
        DdNode* tmp2 = Cudd_bddXnor(manager, int_1[i], int_2[i]);
        Cudd_Ref(tmp2);
        DdNode* new_even = Cudd_bddAnd(manager, tmp2, even);
        Cudd_Ref(new_even);
        Cudd_RecursiveDeref(manager, tmp2);
        Cudd_RecursiveDeref(manager, even);
        even = new_even;
    }
    Cudd_RecursiveDeref(manager, even);
    
    return result;
}

/**Function*************************************************************

  Synopsis    [build a DdNode* type from a function    ]

  Description []

  SideEffects [returned node already referenced]

  SeeAlso     []

***********************************************************************/


DdNode* Simulator::func2node(std::vector<std::string>& func)
{
    assert(!func.empty());
    
    std::stack<DdNode*> waiting;
    for (int i = 0; i < func.size(); ++i)
    {
        if (func[i] == "not")
        {
            assert(waiting.size() >= 1);
            DdNode* temp = waiting.top();
            waiting.pop();
            
            DdNode* new_item = Cudd_Not(temp);
            Cudd_Ref(new_item);
            waiting.push(new_item);
            
            Cudd_RecursiveDeref(manager, temp);
        }
        else if (func[i] == "and" || func[i] == "or" || func[i] == "xor")
        {
            assert(waiting.size() >= 2);
            DdNode* temp_1 = waiting.top();
            waiting.pop();
            DdNode* temp_2 = waiting.top();
            waiting.pop();
             
            DdNode* new_item;              
                 if (func[i] == "and") new_item = Cudd_bddAnd(manager, temp_1, temp_2);
            else if (func[i] == "or")  new_item = Cudd_bddOr (manager, temp_1, temp_2);
            else if (func[i] == "xor") new_item = Cudd_bddXor(manager, temp_1, temp_2);
            Cudd_Ref(new_item);
            waiting.push(new_item);
            
            Cudd_RecursiveDeref(manager, temp_1);
            Cudd_RecursiveDeref(manager, temp_2);
        }
        else
        {
            if (func[i].find_first_not_of( "0123456789" ) == std::string::npos) 
            {
                int ith_var = stoi(func[i]);
                waiting.push(Cudd_bddIthVar(manager, ith_var));
                Cudd_Ref(waiting.top());
            }
            else 
            {
                std::unordered_map<std::string, DdNode*>::iterator iter = defined_var.find(func[i]);
                if (iter == defined_var.end()) 
                {
                    std::cerr << std::endl
                            << "[warning]: Variable \'" << func[i] << "\' is not defined. The variable is taken as constant 0 ..." << std::endl;
                    waiting.push(Cudd_Not(Cudd_ReadOne(manager)));
                    Cudd_Ref(waiting.top());
                }
                else
                {
                    waiting.push(iter->second);
                    Cudd_Ref(waiting.top());
                }
            }
        }
    }
    assert(waiting.size() == 1);
    DdNode *result = waiting.top();
    waiting.pop();
    return result;
}

void Simulator::init_multi_state()
{
    DdNode *tmp, *term, *basis_cube, *label_cube, *entry_cube;

    All_Bdd = new DdNode **[w];
    for (int i = 0; i < w; i++)
        All_Bdd[i] = new DdNode *[r];

    for (int i = 0; i < w; i++)
    {
        for (int j = 0; j < r; j++)
        {
            All_Bdd[i][j] = Cudd_Not(Cudd_ReadOne(manager));
            Cudd_Ref(All_Bdd[i][j]);
        }
    }

    int scale_bits = hasArbitraryInit ? initScaleBits : 0;
    for (int state_index = 0; state_index < multiInitEntries.size(); state_index++)
    {
        label_cube = build_cube(manager, n, multiLabelBits, state_index);
        std::vector<std::string>& entries = multiInitEntries[state_index];
        bool sparse = has_sparse_amplitudes(entries);

        if (entries.size() == 1 && is_bitstring(entries[0]))
        {
            bool ok = true;
            long long real_value = encode_fixed_point(1.0, scale_bits, r, &ok);
            if (!ok)
            {
                std::cerr << "[Error]: Basis initial state for label '" << multiInitLabels[state_index] << "' exceeds the fixed-point range. Increase --r." << std::endl;
                assert(ok);
            }
            basis_cube = build_bitstring_cube(manager, 0, entries[0]);
            entry_cube = Cudd_bddAnd(manager, basis_cube, label_cube);
            Cudd_Ref(entry_cube);
            Cudd_RecursiveDeref(manager, basis_cube);
            for (int bit = 0; bit < r; bit++)
            {
                if (integer_bit(real_value, bit))
                {
                    term = Cudd_bddOr(manager, All_Bdd[3][bit], entry_cube);
                    Cudd_Ref(term);
                    Cudd_RecursiveDeref(manager, All_Bdd[3][bit]);
                    All_Bdd[3][bit] = term;
                }
            }
            Cudd_RecursiveDeref(manager, entry_cube);
        }
        else if (sparse)
        {
            for (std::string& entry : entries)
            {
                size_t sep = entry.find(':');
                if (sep == std::string::npos)
                {
                    std::cerr << "[Error]: Sparse initial-state entry '" << entry << "' for label '" << multiInitLabels[state_index] << "' is missing ':'." << std::endl;
                    assert(sep != std::string::npos);
                }
                std::string basis = entry.substr(0, sep);
                std::string raw_value = entry.substr(sep + 1);
                if (!is_bitstring(basis) || basis.length() != n)
                {
                    std::cerr << "[Error]: Sparse initial-state basis '" << basis << "' for label '" << multiInitLabels[state_index] << "' must be a " << n << "-bit string." << std::endl;
                    assert(is_bitstring(basis) && basis.length() == n);
                }
                bool ok = true;
                std::complex<double> value = parse_complex_value(raw_value, &ok);
                if (!ok)
                {
                    std::cerr << "[Error]: Initial-state amplitude '" << raw_value << "' for label '" << multiInitLabels[state_index] << "' is invalid." << std::endl;
                    assert(ok);
                }
                long long real_value = encode_fixed_point(value.real(), scale_bits, r, &ok);
                long long imag_value = encode_fixed_point(value.imag(), scale_bits, r, &ok);
                if (!ok)
                {
                    std::cerr << "[Error]: Initial-state amplitude for label '" << multiInitLabels[state_index] << "' exceeds the fixed-point range. Increase --r or use smaller amplitudes." << std::endl;
                    assert(ok);
                }
                if (real_value == 0 && imag_value == 0)
                    continue;

                basis_cube = build_bitstring_cube(manager, 0, basis);
                entry_cube = Cudd_bddAnd(manager, basis_cube, label_cube);
                Cudd_Ref(entry_cube);
                Cudd_RecursiveDeref(manager, basis_cube);

                for (int bit = 0; bit < r; bit++)
                {
                    if (integer_bit(real_value, bit))
                    {
                        term = Cudd_bddOr(manager, All_Bdd[3][bit], entry_cube);
                        Cudd_Ref(term);
                        Cudd_RecursiveDeref(manager, All_Bdd[3][bit]);
                        All_Bdd[3][bit] = term;
                    }
                    if (integer_bit(imag_value, bit))
                    {
                        term = Cudd_bddOr(manager, All_Bdd[1][bit], entry_cube);
                        Cudd_Ref(term);
                        Cudd_RecursiveDeref(manager, All_Bdd[1][bit]);
                        All_Bdd[1][bit] = term;
                    }
                }
                Cudd_RecursiveDeref(manager, entry_cube);
            }
        }
        else
        {
            std::vector<std::complex<double>> amplitudes;
            for (std::string& entry : entries)
            {
                bool ok = true;
                std::complex<double> value = parse_complex_value(entry, &ok);
                if (!ok)
                {
                    std::cerr << "[Error]: Initial-state amplitude '" << entry << "' for label '" << multiInitLabels[state_index] << "' is invalid." << std::endl;
                    assert(ok);
                }
                amplitudes.push_back(value);
            }
            for (unsigned long long basis_index = 0; basis_index < amplitudes.size(); basis_index++)
            {
                bool ok = true;
                long long real_value = encode_fixed_point(amplitudes[basis_index].real(), scale_bits, r, &ok);
                long long imag_value = encode_fixed_point(amplitudes[basis_index].imag(), scale_bits, r, &ok);
                if (!ok)
                {
                    std::cerr << "[Error]: Initial-state amplitude for label '" << multiInitLabels[state_index] << "' exceeds the fixed-point range. Increase --r or use smaller amplitudes." << std::endl;
                    assert(ok);
                }
                if (real_value == 0 && imag_value == 0)
                    continue;

                basis_cube = build_cube(manager, 0, n, basis_index);
                entry_cube = Cudd_bddAnd(manager, basis_cube, label_cube);
                Cudd_Ref(entry_cube);
                Cudd_RecursiveDeref(manager, basis_cube);

                for (int bit = 0; bit < r; bit++)
                {
                    if (integer_bit(real_value, bit))
                    {
                        term = Cudd_bddOr(manager, All_Bdd[3][bit], entry_cube);
                        Cudd_Ref(term);
                        Cudd_RecursiveDeref(manager, All_Bdd[3][bit]);
                        All_Bdd[3][bit] = term;
                    }
                    if (integer_bit(imag_value, bit))
                    {
                        term = Cudd_bddOr(manager, All_Bdd[1][bit], entry_cube);
                        Cudd_Ref(term);
                        Cudd_RecursiveDeref(manager, All_Bdd[1][bit]);
                        All_Bdd[1][bit] = term;
                    }
                }
                Cudd_RecursiveDeref(manager, entry_cube);
            }
        }
        Cudd_RecursiveDeref(manager, label_cube);
    }
}

void Simulator::dump_all_bdds(std::string suffix)
{
    if (bddDumpPrefix == "" || manager == nullptr || All_Bdd == nullptr)
        return;

    std::string filename = bddDumpPrefix + "_" + suffix + "_slices.dot";
    FILE *fp = fopen(filename.c_str(), "w");
    if (!fp)
    {
        std::cerr << "[warning]: Cannot open BDD dump file '" << filename << "'." << std::endl;
        return;
    }

    int nVars = Cudd_ReadSize(manager);
    std::vector<std::string> inputNames;
    for (int i = 0; i < nVars; i++)
    {
        if (i < n)
            inputNames.push_back("q" + std::to_string(i));
        else if (i < n + multiLabelBits)
            inputNames.push_back("sel" + std::to_string(i - n));
        else
            inputNames.push_back("aux" + std::to_string(i - n - multiLabelBits));
    }
    std::vector<const char *> inputNamePtrs;
    for (std::string& name : inputNames)
        inputNamePtrs.push_back(name.c_str());

    std::vector<DdNode *> roots;
    std::vector<std::string> outputNames;
    for (int component = 0; component < w; component++)
    {
        for (int bit = 0; bit < r; bit++)
        {
            roots.push_back(All_Bdd[component][bit]);
            outputNames.push_back("c" + std::to_string(component) + "_b" + std::to_string(bit));
        }
    }
    std::vector<const char *> outputNamePtrs;
    for (std::string& name : outputNames)
        outputNamePtrs.push_back(name.c_str());

    Cudd_DumpDot(manager, roots.size(), roots.data(), inputNamePtrs.data(), outputNamePtrs.data(), fp);
    fclose(fp);
    std::cout << "[BDD dump] wrote " << filename << std::endl;
}

void Simulator::dump_big_bdd(std::string suffix)
{
    if (bddDumpPrefix == "" || manager == nullptr || bigBDD == nullptr)
        return;

    std::string filename = bddDumpPrefix + "_" + suffix + "_bigBDD.dot";
    FILE *fp = fopen(filename.c_str(), "w");
    if (!fp)
    {
        std::cerr << "[warning]: Cannot open BDD dump file '" << filename << "'." << std::endl;
        return;
    }

    int nVars = Cudd_ReadSize(manager);
    std::vector<std::string> inputNames;
    for (int i = 0; i < nVars; i++)
    {
        if (i < n)
            inputNames.push_back("q" + std::to_string(i));
        else if (i < n + multiLabelBits)
            inputNames.push_back("sel" + std::to_string(i - n));
        else
            inputNames.push_back("aux" + std::to_string(i - n - multiLabelBits));
    }
    std::vector<const char *> inputNamePtrs;
    for (std::string& name : inputNames)
        inputNamePtrs.push_back(name.c_str());

    DdNode *roots[1] = { bigBDD };
    const char *outputNames[1] = { "bigBDD" };
    Cudd_DumpDot(manager, 1, roots, inputNamePtrs.data(), outputNames, fp);
    fclose(fp);
    std::cout << "[BDD dump] wrote " << filename << std::endl;
}

bool Simulator::has_multi_initial_states()
{
    return !multiInitLabels.empty();
}

std::vector<std::string> Simulator::labels_to_access()
{
    if (selectedInitLabels.empty())
        return multiInitLabels;
    return selectedInitLabels;
}

bool Simulator::select_initial_state(std::string label)
{
    int state_index = -1;
    for (int i = 0; i < multiInitLabels.size(); i++)
    {
        if (multiInitLabels[i] == label)
        {
            state_index = i;
            break;
        }
    }
    if (state_index < 0)
    {
        std::cerr << "[warning]: Initial-state label '" << label << "' is not defined. It is skipped." << std::endl;
        return false;
    }

    baseAllBdd = All_Bdd;
    DdNode *cube = Cudd_ReadOne(manager);
    Cudd_Ref(cube);
    for (int j = multiLabelBits - 1; j >= 0; j--)
    {
        bool bit = (state_index >> j) & 1;
        DdNode *tmp = Cudd_bddAnd(manager, bit ? Cudd_bddIthVar(manager, n + j) : Cudd_Not(Cudd_bddIthVar(manager, n + j)), cube);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, cube);
        cube = tmp;
    }

    All_Bdd = new DdNode **[w];
    for (int i = 0; i < w; i++)
    {
        All_Bdd[i] = new DdNode *[r];
        for (int j = 0; j < r; j++)
        {
            All_Bdd[i][j] = Cudd_Cofactor(manager, baseAllBdd[i][j], cube);
            Cudd_Ref(All_Bdd[i][j]);
        }
    }
    Cudd_RecursiveDeref(manager, cube);
    return true;
}

void Simulator::restore_multi_state()
{
    if (!baseAllBdd)
        return;

    for (int i = 0; i < w; i++)
    {
        for (int j = 0; j < r; j++)
            Cudd_RecursiveDeref(manager, All_Bdd[i][j]);
        delete[] All_Bdd[i];
    }
    delete [] All_Bdd;
    All_Bdd = baseAllBdd;
    baseAllBdd = nullptr;
}

void Simulator::reset_access_state()
{
    if (bigBDD)
    {
        Cudd_RecursiveDeref(manager, bigBDD);
        bigBDD = nullptr;
    }
    for (auto& it: defined_var)
        Cudd_RecursiveDeref(manager, it.second);
    defined_var.clear();
    measure_outcome.clear();
    Node_Table.clear();
    state_count.clear();
    property.clear();
    statevector = "null";
    run_output.clear();
    condition_stack = "";
    range = {0, 0};
    normalize_factor = 1;
}
