# Newtons-Method
Program implementing Newton's method that finds roots of a polynomial of given coefficients

# Usage
Usage: roots -c <coefficients> [OPTIONS]

Options:

   -h --help           Displays this message.

   -V --version        Displays current version.

   -v --verbose        Shows more information about what is being done
                       at the moment.

   -c --coefficients   Sets next numbers as coefficients of your polynomial
                       starting from x^0 and going up, can be floating point.

   -g --guesses        Sets initial guesses from which Newton's method will get
                       more accurate approximations, can be floating point.
                       If none specified, program will try multiple guesses
                       to find all the roots.

   -e --error          Sets largest acceptable error. Default value: 0.0001

Examples:

   roots -c 5 -3 -4 1

   roots -c 5 0 2 1 -g 0 -10 10 30

   roots -v -c 20 5 8 7 1 -e 0.00005
