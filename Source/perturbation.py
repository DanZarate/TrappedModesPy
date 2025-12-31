import numpy as np
import matplotlib.pyplot as plt

def calcular_perfil_fondo(tipo_perturbacion="Rectangular", epsilon=0.5):
    # Parámetros del canal
    lon = 10      # Longitud en x
    ancho = 2    # Ancho en z
    Nx = 500      # Número de puntos en x
    Nz = 100      # Número de puntos en z
    h = 1.0       # Profundidad base del canal

    x = np.linspace(-lon, lon, Nx)
    z = np.linspace(-ancho, ancho, Nz)
    X, Z = np.meshgrid(x, z, indexing='ij')
    
    if tipo_perturbacion == "Rectangular":
        rect_width = 2.0
        H = h * np.ones_like(X)
        mask_band = (np.abs(X) <= rect_width / 2)
        H[mask_band] = epsilon
    elif tipo_perturbacion == "Smooth Crest":
        Am = 1.0
        L0 = 1
        f_x = -Am * (1 / np.cosh(X / L0))**2
        H = h + epsilon * f_x
    elif tipo_perturbacion == "Gaussian":
        V = -np.exp(-X**2)
        H = h + epsilon * V
    else:
        H = h * np.ones_like(X)

    return H


def graficar_comparacion_perturbaciones(epsilon=0.5, mostrar_ecuaciones=True):
    tipos = ["Rectangular", "Smooth Crest", "Gaussian"]
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    colores = ['#1f77b4', '#ff7f0e', '#2ca02c']
    
    for i, tipo in enumerate(tipos):
        H = 1.0 - calcular_perfil_fondo(tipo_perturbacion=tipo, epsilon=epsilon)
        x = np.linspace(-10, 10, H.shape[0])
        ax.plot(x, H[:, H.shape[1]//2], label=f'{tipo}', linewidth=2.5, color=colores[i])
    
    # Línea base
    ax.axhline(y=1.0, color='k', linestyle='--', linewidth=1.5, label='Base Depth $h=1.0$ m')
    
    # Línea de epsilon
    ax.axhline(y=1.0-epsilon, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
    ax.annotate(f'$\\varepsilon = {epsilon}$', 
                xy=(8, 1.0-epsilon), 
                xytext=(8.5, 1.0-epsilon-0.1),
                fontsize=12,
                color='red',
                fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5))
    
    # Título
    ax.set_title(f'Comparison of Channel Depth Profiles ($\\varepsilon = {epsilon}$)', 
                 fontsize=14, fontweight='bold')
    
    # AGREGAR ECUACIONES (versión compatible con matplotlib)
    if mostrar_ecuaciones:
        texto_ecuaciones = (
            r"$\bf{Perturbation\ Equations:}$" + "\n\n"
            r"$\bullet$ Rectangular:" + "\n"
            r"   $H(x) = \varepsilon$ for $|x| \leq 1$" + "\n"
            r"   $H(x) = h$ for $|x| > 1$" + "\n\n"
            r"$\bullet$ Smooth Crest:" + "\n"
            r"   $H(x) = h + \varepsilon A_m \, \mathrm{sech}^2(x/L_0)$" + "\n"
            r"   $(A_m = 1.0,\ L_0 = 1)$" + "\n\n"
            r"$\bullet$ Gaussian:" + "\n"
            r"   $H(x) = h + \varepsilon e^{-x^2}$"
        )
        
        ax.text(0.02, 0.58, texto_ecuaciones,
                transform=ax.transAxes,
                fontsize=9.5,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.85, 
                         edgecolor='navy', linewidth=1.5, pad=0.8))
    
    ax.set_xlabel('$x$ (m)', fontsize=12)
    ax.set_ylabel('Depth $H$ (m)', fontsize=12)
    ax.grid(alpha=0.3, linestyle='--')
    ax.legend(fontsize=10, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(f'trappedmodes\\Output\\channel_depth_profiles_comparison_eps{epsilon}.png', 
                dpi=300, bbox_inches='tight')
    print(f'Figura guardada: channel_depth_profiles_comparison_eps{epsilon}.png')


# Probar con diferentes valores de epsilon
graficar_comparacion_perturbaciones(epsilon=0.5, mostrar_ecuaciones=True)

# Si quieres generar múltiples gráficas:
# for eps in [0.3, 0.5, 0.7]:
#     graficar_comparacion_perturbaciones(epsilon=eps, mostrar_ecuaciones=True)