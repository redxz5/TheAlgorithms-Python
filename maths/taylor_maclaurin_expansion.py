"""
Taylor Series : https://en.wikipedia.org/wiki/Taylor_series
Maclaurin Series : Taylor Series at point x = 0.
"""


from sympy import symbols, factorial, diff, Add, Symbol
"""
Variable "var" can be declared using symbols("var")
Functions like "sin" "log" "exp" need to be imported to be used
For display of Expansion terms import "pprint"
"""
x = symbols('x')
y = symbols('y')


def taylor_expansion(func: Add, point: int, terms: int, var: Symbol):
    """
    Returns Taylor expansion of the "func" up to the required number of terms.
    Parameters :
    func  - The function to be approximated
    point - Point at which the value is to be approximated
    terms - Number of terms to consider (more terms increases accuracy)
    var   - The independent variable with respect to which
            the Taylor expansion is calculated.

    >>> from sympy import pprint, exp, log, sin
    >>> pprint(taylor_expansion(exp(x), 0, 5, x))
     4    3    2        
    x    x    x         
    ── + ── + ── + x + 1
    24   6    2         

    >>> pprint(taylor_expansion(log(x + 1), 0, 5, x))
     5    4    3    2    
    x    x    x    x     
    ── - ── + ── - ── + x
    5    4    3    2     

    >>> pprint(taylor_expansion(sin(x), 0, 5, x))
     5     4    3    2    
    x     x    x    x     
    ─── - ── + ── - ── + x
    120   24   6    2     

    >>> pprint(taylor_expansion(exp(x * y), 1, 5, x))
     4        4  y    3        3  y    2        2  y                    
    y ⋅(x - 1) ⋅ℯ    y ⋅(x - 1) ⋅ℯ    y ⋅(x - 1) ⋅ℯ               y    y
    ────────────── + ────────────── + ────────────── + y⋅(x - 1)⋅ℯ  + ℯ 
          24               6                2                           

    >>> pprint(taylor_expansion(exp(x * y), 1, -1, x))
    3rd argument is for no. of terms, provide a natural number
    <BLANKLINE>
    """
    t_y = symbols('t_y')
    expansion = func.subs(var, point)
    d = func
    if (expansion == 0):
        terms += 1
    try:
        k = 1
        while (k < terms):
            d = diff(d, var)
            term = (d * ((t_y - point) ** k)) / factorial(k)
            term = term.subs(var, point)
            if (term == 0):
                continue
            term = term.subs(t_y, var)
            expansion = term + expansion
            k += 1
            if (d == 0 and k < terms):
                print("only ", k - 1, " terms present")
        if (terms < 1):
            print("3rd argument is for no. of terms, provide a natural number")
            return ''
        expansion = expansion.subs(t_y, var)
        return expansion
    except TypeError:
        print("3rd argument denotes number of terms, provide a natural number")
        return ''


def taylor_value(func: Add, point: int, terms: int, var: Symbol):
    """
    Substitutes "var" in order to calcualate value of taylor approximation.
    Parameters :
    func  - The function to be approximated
    point - Point at which the value is to be approximated
    terms - Number of terms to consider (more terms increases accuracy)
    var   - The independent variable with respect to which
            the Taylor expansion is calculated.

    >>> from sympy import pprint, exp, log, sin
    >>> print(taylor_value(exp(x), 1, 10, x))
    2.71828182845905
    """
    f = taylor_expansion(func, point, terms, var)
    return f.evalf(subs={var: point})


def maclaurin_expansion(func: Add, terms: int, var: Symbol):
    """
    Returns maclaurin expansion of "func" up to the required number of terms.
    Parameters :
    func  - The function to be approximated
    terms - Number of terms to consider (more terms increases accuracy)
    var   - The independent variable with respect to which
            the Taylor expansion is calculated.

    >>> from sympy import pprint, exp, log, sin
    >>> pprint(maclaurin_expansion(exp(x), 5, x))
     4    3    2        
    x    x    x         
    ── + ── + ── + x + 1
    24   6    2         

    >>> pprint(maclaurin_expansion(log(x + 1), 5, x))
     5    4    3    2    
    x    x    x    x     
    ── - ── + ── - ── + x
    5    4    3    2     

    >>> pprint(maclaurin_expansion(sin(x), 5, x))
     5     4    3    2    
    x     x    x    x     
    ─── - ── + ── - ── + x
    120   24   6    2     

    """
    return taylor_expansion(func, 0, terms, var)


def maclaurin_value(func: Add, terms: int, var: Symbol):
    """
    Substitutes "var" in order to calcualate value of maclaurin approximation.
    Parameters :
    func  - The function to be approximated
    terms - Number of terms to consider (more terms increases accuracy)
    var   - The independent variable with respect to which
            the Taylor expansion is calculated.

    >>> from sympy import pprint, exp, log, sin
    >>> print(maclaurin_value(exp(1+x), 5, x))
    2.71828182845905
    """
    return taylor_value(func, 0, terms, var)
