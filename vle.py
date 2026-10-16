import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

def antoine_inverse(P, A, B, C):
    """
    Εύρεση θερμοκρασίας βρασμου (σε °C) για μια δεδομένη τάση ατμού P (σε mmHg) χρησιμοποιώντας την εξίσωση Antoine.
    """
    T = (B / (A - np.log10(P))) - C
    return T

def antoine(T, A, B, C):
    """
    Εύρεση πίεσης ατμού (σε mmHg) για μια δεδομένη θερμοκρασία T (σε °C) χρησιμοποιώντας την εξίσωση Antoine.
    """
    P = 10 ** (A - (B / (T + C)))
    return P

def bubble_residual(T, x1, P, comp1, comp2, vl=None):
    """
    Υπολογισμός του υπολοίπου της εξίσωσης για το σημείο βρασμού ενός δυαδικού μίγματος.
    """

    p1 = antoine(T, *comp1)
    p2 = antoine(T, *comp2)

    if vl is None:
        g1 = g2 = 1.0
    else:
        g1, g2 = van_laar(x1, *vl)

    return x1*g1*p1 + (1-x1)*g2*p2 - P

def bubble_point_temperature(x1, P, comp1, comp2, vl=None):

    """
    Υπολογισμός της θερμοκρασίας βρασμού ενός δυαδικού μίγματος μεσω Antoine.
    """
    T_low = antoine_inverse(P, *comp1)  # Κατω ακρο
    T_high = antoine_inverse(P, *comp2)  # Ανω ακρο
    T_bubble = brentq(bubble_residual, T_low-20, T_high+20, args=(x1, P, comp1, comp2, vl))
    return T_bubble

def bubble_point(x1, P, comp1, comp2, vl=None):

    """
    Υπολογισμός της θερμοκρασίας βρασμού ενός δυαδικού μίγματος μεσω Antoine.
    """
    T_bubble = bubble_point_temperature(x1, P, comp1, comp2, vl )
    if vl is None:
        g1 = 1.0
    else:
        g1, _ = van_laar(x1, *vl)
    y1 = (x1 * g1 * antoine(T_bubble, *comp1)) / P
    return T_bubble, y1

def van_laar(x1, A12, A21):
    

    if x1 < 0.0 or x1 > 1:
        raise ValueError("x1 must be between 0 and 1.")
    x2 = 1 - x1
    if x1 == 0:
        return np.exp(A12), 1.0
    if x2 == 0:
        return 1.0, np.exp(A21)
    

    d1 = 1.0 + (A12 * x1) / (A21 * x2)
    d2 = 1.0 + (A21 * x2) / (A12 * x1)
    
    return np.exp(A12/d1**2), np.exp(A21/d2**2)

    




def main():
   water = (8.07131, 1730.63, 233.426)
   methanol = (8.08097, 1582.271, 239.726)
   VL = (0.8407, 0.5386)

   x1 = float(input("Enter the mole fraction of methanol in the liquid phase (x1): "))
   P = float(input("Enter the total pressure in mmHg (P): "))

   if not 400 <= P <= 1500:
       print("Warning: Pressure is outside the typical range for this system (400-1500 mmHg).")

   T_id, y_id = bubble_point(x1, P, methanol, water)
   T_vl, y_vl = bubble_point(x1, P, methanol, water, vl=VL)

   print()
   print(f"At P = {P:.2f} mmHg and x1 = {x1:.3f}:")
   print(f" Ideal (Raoult) : T = {T_id:7.2f} °C, y1 = {y_id:.4f}")
   print(f" van Laar       : T = {T_vl:7.2f} °C, y1 = {y_vl:.4f}")
   print(f" Difference    : ΔT = {T_vl - T_id:+6.2f} °C, Δy1 = {y_vl - y_id:.4f}")


   x_vals = np.linspace(0, 1, 101)
   y_ideal, y_real = [], []
   T_ideal, T_real = [], []
   
   for x in x_vals:
      Ti, yi = bubble_point(x, P, methanol, water)
      Tr, yr = bubble_point(x, P, methanol, water, vl=VL)
      T_real.append(Tr)
      T_ideal.append(Ti)
      y_real.append(yr)
      y_ideal.append(yi)

   plt.figure(figsize=(6, 6))
   plt.plot(x_vals, y_ideal, '--', label='Ideal (Raoult)', color='gray')
   plt.plot(x_vals, y_real, '-', label='van Laar', color='blue')
   plt.plot([0, 1], [0, 1], "k:", label = "y = x")
   plt.xlabel('Liquid Mole Fraction of Methanol (x1)')
   plt.ylabel('Vapor Mole Fraction of Methanol (y1)')
   plt.title(f'Vapor-Liquid Equilibrium Curve for Methanol-Water System at {P:.0f} mmHg')
   plt.xlim(0, 1)
   plt.ylim(0, 1)
   plt.legend()
   plt.grid(alpha=0.3)
   plt.savefig('Equilibrium_Curve_.png', dpi=150, bbox_inches='tight')
   plt.show()

   plt.figure(figsize=(7, 5))
   plt.plot(x_vals, T_ideal, '--', color='gray', label='Ideal (Raoult)')
   plt.plot(x_vals, T_real, '-', color='blue', label='van Laar')
   plt.plot(y_ideal, T_ideal, '--', color='gray')
   plt.plot(y_real, T_real, '-', color='blue')
   plt.xlabel('Liquid Mole Fraction of Methanol (x1/y1)')
   plt.ylabel('Temperature (°C)')
   plt.title(f'Bubble Point Temperature for Methanol-Water System at {P:.0f} mmHg')
   plt.legend()
   plt.grid(alpha=0.3)
   plt.savefig('Bubble_Point_Temperature_.png', dpi=150, bbox_inches='tight')
   plt.show()
   
if __name__ == "__main__":
    main()