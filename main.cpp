#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <algorithm>

#include <cmath>

const std::string version = "Alpha v1.0";

struct Data
{
    char** arguments;
    std::vector<double> coefficients;   //from x^0 to the highest power
    std::vector<double> guesses;

    std::unique_ptr<std::vector<double>> roots;

    double acceptable_error = 0.0001;
    
    int max_iterations = 0xFFFF;

    bool verbose = false;

};

Data data;

//processing functions
void process_input(int argc, char* argv[]);
std::string vec_str(const std::vector<double> &vec);

//numerical/analytical functions
std::vector<double>* derivative(const std::vector<double> &f);
std::vector<double>* get_roots(const std::vector<double> &f);
std::vector<double>* newton_method(const std::vector<double> &f, const std::vector<double> &d, const std::vector<double> &guesses);
double evaluate(const std::vector<double> &f, double x);

//argument functions
void display_help(int &index);
void get_coefficients(int &index);
void get_guesses(int &index);
void display_version(int &index);
void set_verbose(int &index);
void set_error(int &index);

int main(int argc, char* argv[])
{
    process_input(argc, argv);

    if(data.guesses.empty()) data.roots = std::unique_ptr<std::vector<double>>(get_roots(data.coefficients));
    else
    {
        auto &der = *derivative(data.coefficients);
        data.roots = std::unique_ptr<std::vector<double>>(newton_method(data.coefficients, der, data.guesses));
        delete &der;
    }

    if(data.roots->empty())
    {
        std::cout << "Sorry! No real roots found!\n";
        return 0;
    }

    for(int i = 0; i < data.roots->size(); i++)
    {   
        std::cout << 'x' << i << " = " << data.roots->at(i) << '\n';
    }

    return 0;
}


void process_input(int argc, char* argv[])
{
    data.arguments = argv;

    std::unordered_map<std::string, void (*)(int&)> arg_map = {
        {"-h", display_help},       {"--help", display_help},       
        {"-c", get_coefficients},   {"--coefficients", get_coefficients},
        {"-g", get_guesses},        {"--guessees", get_guesses},    
        {"-V", display_version},    {"--version", display_version},     
        {"-v", set_verbose},        {"--verbose", set_verbose},
        {"-e", set_error},          {"--error", set_error}
    };


    //note that index may change, as some functions go through next arguments
    int dummy = 0;
    int index = 1;

    if(argc == 1) display_help(dummy);

    while(index < argc)
    {   
        try
        {
            //try to call function from the map
            arg_map.at(argv[index])(index);
        }
        catch(const std::exception& e)
        {
            //invalid argument
            display_help(index);
        }
        index++;
    }

    if(data.coefficients.empty()) std::exit(0);

    //get rid of 0s near the highest powers
    while(!*(data.coefficients.end() - 1))
    {
        data.coefficients.pop_back();
    }

    //display polynomial in human-friendly format
    std::cout << "Your polynomial: ";

    for(int i = data.coefficients.size() - 1; i >= 0; i--)
    {
        if(data.coefficients.at(i) == 0) continue;

        if(data.coefficients.at(i) < 0) std::cout << " - ";
        //if it's not the first number
        else if(i < data.coefficients.size() - 1)  std::cout << " + ";
        
        if(std::fabs(data.coefficients.at(i) != 1) || i == 0) std::cout << std::fabs(data.coefficients.at(i)) << ' ';
        
        if(i > 0) std::cout << 'x';

        if(i > 1) std::cout << '^' << i;
    }
    std::cout << "\n\n";

    if(!data.verbose) return;
    
    std::cout << "Coefficients: { ";
    for(double x : data.coefficients) std::cout << x <<',';
    std::cout << "\b }\nGuesses: { ";
    for(double x : data.guesses) std::cout << x << ',';
    std::cout << "\b }\n\n";

    std::cout << "Acceptable error: " << data.acceptable_error << "\n\n";
    
}

std::vector<double>* derivative(const std::vector<double> &f)
{
    auto *d = new std::vector<double>;

    for(int i = 1; i < f.size(); i++)
    {
        d->push_back(i * f.at(i));
    }

    return d;
}

double evaluate(const std::vector<double> &f, double x)
{
    double result = 0;

    for(int i = 0; i < f.size(); i++)
        result += std::pow(x, i) * f.at(i);

    return result;
}

std::vector<double>* get_roots(const std::vector<double> &f)
{
    if(data.verbose) std::cout << "Getting roots of " << vec_str(f) << "\n";

    auto *result = new std::vector<double>;

    if(f.size() < 2) return result; 

    if(f.size() == 2)
    {
        result->push_back( - f.at(0) / f.at(1));
        return result;
    }

    if(f.size() == 3)
    {
        double delta = f.at(1) * f.at(1) - 4 * f.at(2) * f.at(0);
        
        if(delta < 0) return result;

        double sqrt = std::sqrt(delta);

        result->push_back( (- f.at(1) - sqrt) / (2 * f.at(2)));

        if(delta) result->push_back( (- f.at(1) + sqrt) / (2 * f.at(2)));

        return result;
    }

    //get guesses
    std::unique_ptr<std::vector<double>> dfodx(derivative(f));
    std::unique_ptr<std::vector<double>> extrema(get_roots(*dfodx));    //RECURSION ALERT

    if(data.verbose)
    {
        std::cout << "Continuing work with " << vec_str(f) << '\n' 
        << "Extrema: " << vec_str(*extrema) << '\n';
    }

    if(extrema->empty())
        return newton_method(f, *dfodx, std::vector<double>{0});

    if(extrema->size() == 1)
        return newton_method(f, *dfodx, std::vector<double>{extrema->at(0) - 1, extrema->at(0) + 1});

    std::vector<double> guesses;

    //before first
    guesses.push_back(0.5 * (3 * extrema->at(0) - extrema->at(1)));

    //in between
    for(int i = 0; i < extrema->size() - 1; i++ )
        guesses.push_back(0.5 * (extrema->at(i) + extrema->at(i + 1)));
        
    //after the last
    guesses.push_back(0.5 * (3 * extrema->at(extrema->size() - 1) - extrema->at(extrema->size() - 2)));

    return newton_method(f, *dfodx, guesses);
}

std::vector<double>* newton_method(const std::vector<double> &f, const std::vector<double> &d, const std::vector<double> &guesses)
{
    auto *roots = new std::vector<double>;

    for(auto guess : guesses)
    {
        double lastx = guess;
        double y = evaluate(f, lastx);
        double x = lastx - y / evaluate(d, lastx);

        for(int i = 0; y != 0 && i < data.max_iterations; i++ )
        {
            lastx = x;
            y = evaluate(f, x);
            x -= y / evaluate(d, x);
        }

        if(std::find(roots->begin(), roots->end(), x) == roots->end() && std::fabs(y) < data.acceptable_error)
            roots->push_back(x);
    }

    return roots;
}

void get_coefficients(int &index)
{
    double c;
    while(true)
    {
        try
        {
            c = std::stof(data.arguments[++index]);
        }
        catch(const std::exception& e)
        { 
            index--;
            break; 
        }
        
        data.coefficients.push_back(c);
    }
}

void get_guesses(int &index)
{
    double c;
    while(true)
    {
        try
        {
            c = std::stof(data.arguments[++index]);
        }
        catch(const std::exception& e)
        { 
            index--;
            break; 
        }
        
        data.guesses.push_back(c);
    }
}

void display_help(int &index)
{
    if(data.arguments[index] == std::string("-h") || data.arguments[index] == std::string("--help"))
        std::cout << "It is a program that uses Newton's method to approximate roots of a polynomial.\n";
    else if(index)
        std::cout << "Unrecognized option (index: " << index << ")!\n\n";

    std::cout << "Usage: roots -c <coefficients> [OPTIONS]\n\n"

    "Options:\n\n"

    "   -h --help           Displays this message.\n\n"

    "   -V --version        Displays current version.\n\n"

    "   -v --verbose        Shows more information about what is being done\n"
    "                       at the moment.\n\n"

    "   -c --coefficients   Sets next numbers as coefficients of your polynomial\n"
    "                       starting from x^0 and going up, can be floating point.\n\n"

    "   -g --guesses        Sets initial guesses from which Newton's method will get\n"
    "                       more accurate approximations, can be floating point.\n"
    "                       If none specified, program will try multiple guesses\n"
    "                       to find all the roots.\n\n"

    "   -e --error          Sets largest acceptable error. Default value: 0.0001\n\n"
    
    "Examples:\n\n"

    "   roots -c 5 -3 -4 1\n\n"

    "   roots -c 5 0 2 1 -g 0 -10 10 30\n\n"
    
    "   roots -v -c 20 5 8 7 1 -e 0.00005\n";
}

void display_version(int &index)
{
    std::cout << "Newton's method program made by Maciej Suski\n"
    "License: THE BEER-WARE LICENSE. It means you can use it, but if we meet,\n"
    "         you buy me a beer.\n"
    "Version: " << version << std::endl;
}

void set_verbose(int &index)
{ data.verbose = true; }

void set_error(int &index)
{
    try
    {
        data.acceptable_error = std::stod(data.arguments[++index]);
    }
    catch(const std::exception& e)
    {
        std::cout << "Invalid argument (index: " << index << "), expected number.\n";
        std::exit(-1);
    }
}

std::string vec_str(const std::vector<double> &vec)
{
    std::string str = "{ ";

    if(vec.empty())
        return str + "}";

    for(int i = 0; i < vec.size() - 1; i++)
        str += std::to_string(vec.at(i)) + ", ";
 
    str += std::to_string(vec.at(vec.size() - 1));

    return str + " }";
}