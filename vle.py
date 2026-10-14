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

def bubble_residual(T, x1, P, comp1, comp2):
    """
    Υπολογισμός του υπολοίπου της εξίσωσης για το σημείο βρασμού ενός δυαδικού μίγματος.
    """

    p1 = antoine(T, *comp1)
    p2 = antoine(T, *comp2)

    return x1*p1 + (1-x1)*p2 - P 

def bubble_point_temperature(x1, P, comp1, comp2):
    """
    Υπολογισμός της θερμοκρασίας βρασμού ενός δυαδικού μίγματος μεσω Antoine.
    """
    T_low = antoine_inverse(P, *comp1)  # Κατω ακρο
    T_high = antoine_inverse(P, *comp2)  # Ανω ακρο
    T_bubble = brentq(bubble_residual, T_low-1, T_high+1, args=(x1, P, comp1, comp2))
    return T_bubble

def bubble_point(x1, P, comp1, comp2):
    """
    Υπολογισμός της θερμοκρασίας βρασμού ενός δυαδικού μίγματος μεσω Antoine.
    """
    T_bubble = bubble_point_temperature(x1, P, comp1, comp2)
    y1 = (x1 * antoine(T_bubble, *comp1)) / P
    return T_bubble, y1

def main():
   water = (8.07131, 1730.63, 233.426)
   methanol = (8.08097, 1582.271, 239.726)
   print(bubble_point_temperature(0.5, 760, methanol, water))
   T, y1 = bubble_point(0.5, 760, methanol, water)
   print(f"Bubble point temperature: {T:.2f} °C, Vapor mole fraction of component 1: {y1:.4f}")
   x_vals = np.linspace(0, 1, 11)
   y_vals = []
   T_vals = []
   for x1 in x_vals:
      T, y1 = bubble_point(x1, 760, methanol, water)
      T_vals.append(T)
      y_vals.append(y1)
   plt.figure(figsize=(6, 6))
   plt.plot(x_vals, y_vals, label='Bubble Point Temperature (°C)', color='blue')
   plt.plot([0, 1], [0, 1], "k--", label = "y=x")
   plt.xlabel('Liquid Mole Fraction of Methanol (x1)')
   plt.ylabel('Vapor Mole Fraction of Methanol (y1)')
   plt.title('Vapor-Liquid Equilibrium Curve for Methanol-Water System at 760 mmHg')
   plt.legend(["Equilibrium Curve", "y=x Line"])
   plt.grid(alpha=0.3)
   plt.show()
    

main()