# Python3 program to find modulo  
# of floating point numbers. 
  
def findMod(a, b): 
  
    # Handling negative values 
    if (a < 0): 
        a = -a 
    if (b < 0): 
        b = -b 
  
    # Finding mod by repeated subtraction 
    mod = a 
    while (mod >= b): 
        mod = mod - b 
  
    # Sign of result typically  
    # depends on sign of a. 
    if (a < 0): 
        return -mod 
  
    return mod 
  
# Driver code 
if __name__ == "__main__":
    a = 9.7; b = 2.3
    print(findMod(a, b)) 
  
# This code is contributed by Anant Agarwal. 
