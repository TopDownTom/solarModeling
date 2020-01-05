# Runge-Kutta 4 
# Uncomment any of the print statements to see the output at each step
def rk4(function,y,r,step):
    F1 = function(y,r)
#    print("F1: {}".format(F1))
    F2 = function( y + step*F1/2 , r + .5*step)
#    print("F2: {}".format(F2))
    F3 = function( y + step*F2/2 , r + .5*step)
#    print("F3: {}".format(F3))
    F4 = function( y + step*F3 , r + step)
#    print("F4: {}".format(F4))
    FN = y + ( step*( F1 + F4 + 2*(F2 + F3) )/6 )
#    print("FN: {}".format(FN))
#    print("\n______")
    return FN
