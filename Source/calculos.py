# calculos.py
import numpy as np
condicionInicial = None
condicionFrontera = None   
H = None
zoom = False
X = None
Z = None
x = None
control = None
ganancia = 0.0
coordenadas_x_z = (0.0, 0.0)

def calcular_canal(lon, ancho, Nx=100, Nz=100, h=1.0, T=100.0, condicion_inicial=None, condicion_frontera=None, tipo_perturbacion=None, control_tipo=None, ganancia_val=0.0, coordenadas=(0.0, 0.0)):
    global condicionInicial
    global condicionFrontera
    global H
    global zoom
    global X, Z, x
    global control
    global ganancia
    global coordenadas_x_z
    coordenadas_x_z = coordenadas
    control = control_tipo
    ganancia = ganancia_val

    condicionInicial = condicion_inicial
    condicionFrontera = condicion_frontera
    
    # Si alguna condicion de frontera es 'Infinite', activar zoom
    if condicion_frontera is not None:
        if 'Infinite' in condicion_frontera:
            zoom = True
        else:
            zoom = False

    print(f"Total Simulation Time T: {T}")

    # Gravedad estándar
    g = 9.81  # m/s^2

    # Cálculo de los pasos espaciales
    dx = 2*lon / (Nx - 1)
    dz = 2*ancho / (Nz - 1)
    print(f"Spatial step in x (dx): {dx}")
    print(f"Spatial step in z (dz): {dz}")

    dt = 0.9 * min(dx, dz) / np.sqrt(2 * g * h)
    print(f"Time step (dt): {dt}")
    Nt = int(T / dt) + 1
    print(f"Number of time steps (Nt): {Nt}")

    x = np.linspace(-lon, lon, Nx)
    z = np.linspace(-ancho, ancho, Nz)
    X, Z = np.meshgrid(x, z, indexing='ij')
    
    epsilon = 0.5

    if tipo_perturbacion == "Rectangular":
        rect_width = 2.0   # Ancho en dirección x (aumentado para mejor visualización)
        H = h * np.ones_like(X)
        mask_band = (np.abs(X) <= rect_width / 2)
        H[mask_band] = epsilon # Profundidad reducida en la banda central
    elif tipo_perturbacion == "Smooth Crest":
        Am = 1.0         # Amplitud de la perturbación
        L0 = 1.0       # Longitud característica (ancho de la cresta)
        f_x = -Am * (1 / np.cosh(X / L0))**2
        H = h + epsilon * f_x
    elif tipo_perturbacion == "Gaussian":
        V = -np.exp(-X**2)
        H = h + epsilon * V
    else:
        H = h * np.ones_like(X)

    frames, E_history, magnitude_spectrum, kx_freqs, kz_freqs, amplitud_en_coordenada = simular_evolucion(
        dx, dz, dt, Nx, Nz, Nt, g)
    
    return frames.copy(), X, Z, E_history, dt, magnitude_spectrum, kx_freqs, kz_freqs, amplitud_en_coordenada

def compute_laplacian(psi, dx, dz):
    global H
    laplacian = np.zeros_like(psi)
    # Promedios de H en las caras (x-dirección)
    Hx_plus = 0.5 * (H[2:, 1:-1] + H[1:-1, 1:-1])
    Hx_minus = 0.5 * (H[1:-1, 1:-1] + H[:-2, 1:-1])
    # Término en x usando forma conservativa
    laplacian[1:-1, 1:-1] += (Hx_plus * (psi[2:, 1:-1] - psi[1:-1, 1:-1]) - 
                              Hx_minus * (psi[1:-1, 1:-1] - psi[:-2, 1:-1])) / dx**2
    # Promedios de H en las caras (z-dirección)
    Hz_plus = 0.5 * (H[1:-1, 2:] + H[1:-1, 1:-1])
    Hz_minus = 0.5 * (H[1:-1, 1:-1] + H[1:-1, :-2])
    # Término en z usando forma conservativa
    laplacian[1:-1, 1:-1] += (Hz_plus * (psi[1:-1, 2:] - psi[1:-1, 1:-1]) - 
                              Hz_minus * (psi[1:-1, 1:-1] - psi[1:-1, :-2])) / dz**2
    return laplacian

def condiciones_iniciales(X, Z):
    global condicionInicial
    psi0 = np.zeros_like(X)
    if condicionInicial == 'None':
        psi0 = np.zeros_like(X)
        return psi0
    elif condicionInicial == 'Water Drop':
        psi0 = np.exp(-(X)**2) * np.cos(np.pi * Z) * 0.001
        return psi0
    elif condicionInicial == 'Sinusoidal':
        kx = np.pi / (X.max() - X.min())
        kz = np.pi / (Z.max() - Z.min())
        psi0 = np.sin(kx * X) * np.sin(kz * Z)
        return psi0
    elif condicionInicial == 'Gaussian':
        psi0 = np.exp(-((X/10)**2 + (Z/5)**2))
        return psi0
    else:
        print("Condición inicial no reconocida. Usando Gaussian por defecto.")
        # Ejemplo de condición inicial: una perturbación Gaussian
        psi0 = np.exp(-((X/10)**2 + (Z/5)**2))
        return psi0
    return psi0

# Aplicar condiciones de frontera Dirichlet, Neumann, absorting, Open, etc.
def aplicar_condiciones_frontera(psi_np1, psi_n, c, dt, dx):
    global condicionFrontera
    frontera_lat_izq = condicionFrontera[0]  # Lateral izq.
    frontera_lat_der = condicionFrontera[1]  # Lateral der.
    frontera_lat_sup = condicionFrontera[2]  # Lateral sup.
    frontera_lat_inf = condicionFrontera[3]  # Lateral inf.

    # x = -L
    if frontera_lat_izq == 'Dirichlet':
        psi_np1[0, :] = 0
    elif frontera_lat_izq == 'Neumann':
        psi_np1[0, :] = psi_np1[1, :]
    elif frontera_lat_izq == 'Open' or frontera_lat_izq == 'Infinite':
        psi_np1[0, :] = psi_n[1, :] + (c * dt - dx) / (c * dt + dx) * (psi_np1[-2, :] - psi_n[-1, :])
        psi_np1[0, :] *= 0.8
        psi_np1[1, :] *= 0.9
    elif frontera_lat_izq == 'Absorbing':
        psi_np1[0, :] = 0.9 * psi_np1[1, :]  
    elif frontera_lat_izq == 'Absorbing2':
        psi_np1[0, :] = 0.8 * psi_np1[1, :]

    # x = +L
    if frontera_lat_der == 'Dirichlet':
        psi_np1[-1, :] = 0
    elif frontera_lat_der == 'Neumann':
        psi_np1[-1, :] = psi_np1[-2, :]
    elif frontera_lat_der == 'Open' or frontera_lat_der == 'Infinite':
        psi_np1[-1, :] = psi_n[-2, :] + (c * dt - dx) / (c * dt + dx) * (psi_np1[1, :] - psi_n[1, :])
        psi_np1[-1, :] *= 0.8
        psi_np1[-2, :] *= 0.9
    elif frontera_lat_der == 'Absorbing':
        psi_np1[-1, :] = 0.9 * psi_np1[-2, :]  
    elif frontera_lat_der == 'Absorbing2':
        psi_np1[-1, :] = 0.8 * psi_np1[-2, :]

    # z = -d
    if frontera_lat_inf == 'Dirichlet':  
        psi_np1[:, 0] = 0
    elif frontera_lat_inf == 'Neumann':
        psi_np1[:, 0] = psi_np1[:, 1]
    elif frontera_lat_inf == 'Open':
        psi_np1[:, 0] = psi_n[:, 1] + (c * dt - dx) / (c * dt + dx) * (psi_np1[:, -2] - psi_n[:, -1])
        psi_np1[:, 0] *= 0.8
        psi_np1[:, 1] *= 0.9
    elif frontera_lat_inf == 'Absorbing':
        psi_np1[:, 0] = 0.9 * psi_np1[:, 1]
    elif frontera_lat_inf == 'Absorbing2':
        psi_np1[:, 0] = 0.8 * psi_np1[:, 1]
    
    # z = +d
    if frontera_lat_sup == 'Dirichlet':
        psi_np1[:, -1] = 0
    elif frontera_lat_sup == 'Neumann':
        psi_np1[:, -1] = psi_np1[:, -2]
    elif frontera_lat_sup == 'Open':
        psi_np1[:, -1] = psi_n[:, -2] + (c * dt - dx) / (c * dt + dx) * (psi_np1[:, 1] - psi_n[:, 1])
        psi_np1[:, -1] *= 0.8
        psi_np1[:, -2] *= 0.9
    elif frontera_lat_sup == 'Absorbing':
        psi_np1[:, -1] = 0.9 * psi_np1[:, -2]
    elif frontera_lat_sup == 'Absorbing2':
        psi_np1[:, -1] = 0.8 * psi_np1[:, -2]

    return psi_np1

def actualizar_psi(psi, psi_prev, laplacian, dt, g):
    global control
    global ganancia

    # Actualizar psi usando un esquema explícito
    psi_new = (2 * psi - psi_prev + dt**2 * g * laplacian)
    
    if control == "Velocity Feedback":
        velocidad = (psi - psi_prev) / dt
        psi_new -= ganancia * velocidad * dt**2 

    return psi_new

def simular_evolucion(dx, dz, dt, Nx, Nz, Nt, g):

    global H
    global X, Z, x
    global zoom
    global coordenadas_x_z
    #Print media de H
    print(f"Profundidad media H: {np.mean(H)}")

    frames = []

    # Inicializar Energía
    E_history = []
    E_cin_history = []
    E_pot_history = []
    frames_transformada = []
    kx_freqs = []  # Para almacenar las frecuencias
    kz_freqs = []

    # Inicializar Amplitud de la onda en la coordenada dada
    x0, z0 = coordenadas_x_z
    idx_x = (np.abs(x - x0)).argmin()
    idx_z = (np.abs(Z[0, :] - z0)).argmin()
    amplitud_en_coordenada = []

    psi_nm1 = np.zeros((Nx, Nz))  # psi en el tiempo n-1
    psi_n = np.zeros((Nx, Nz))    # psi en el tiempo n
    psi_np1 = np.zeros((Nx, Nz))  # psi en el tiempo n+1

    psi_n = condiciones_iniciales(X, Z)
    psi_nm1 = psi_n.copy()  # Asumimos que la condición inicial en n-1 es igual a n

    # Determinar la región para la transformada
    if zoom:
        xmin, xmax = -10, 10
        zona_x = (x >= xmin) & (x <= xmax)
        zona_idx = np.where(zona_x)[0]
        Nx_fft = len(zona_idx)
    else:
        zona_idx = None
        Nx_fft = Nx

    # Calcular las frecuencias espaciales (una sola vez)
    kx = np.fft.fftshift(np.fft.fftfreq(Nx_fft, dx))  # ciclos por unidad de longitud
    kz = np.fft.fftshift(np.fft.fftfreq(Nz, dz))

    for n in range(1, Nt):
        laplacian = compute_laplacian(psi_n, dx, dz)
        psi_np1 = actualizar_psi(psi_n, psi_nm1, laplacian, dt, g)
        psi_np1 = aplicar_condiciones_frontera(psi_np1, psi_n, np.sqrt(g * H.mean()), dt, dx)

        # ---- Energía ----
        E_total, E_cin, E_pot, energy_density = calcular_energia_total(
            psi_n, psi_nm1, dx, dz, dt, g
        )
        
        E_history.append(E_total)
        E_cin_history.append(np.sum(E_cin) * dx * dz)
        E_pot_history.append(np.sum(E_pot) * dx * dz)
        #-- Fin energía ----

        # Amplitud en la coordenada dada
        amplitud_en_coordenada.append(psi_n[idx_x, idx_z])

        # --- Transformada de Fourier 2D ---
        # Aplicar en la región de zoom si está activado
        if zoom and zona_idx is not None and len(zona_idx) > 0:
            psi_fft = psi_n[zona_idx, :]
        else:
            psi_fft = psi_n
        
        fft_psi = np.fft.fft2(psi_fft)
        fft_shifted = np.fft.fftshift(fft_psi)
        magnitude_spectrum = np.abs(fft_shifted)

        # Espectro de potencia (mejor que magnitud simple)
        power_spectrum = np.abs(fft_shifted)**2
        # ----------------------------------

        # Almacenar algunos frames para visualización
        if n % 10 == 0: #(Nt // 10) == 0:
            frames.append(psi_n.copy())
            #frames_transformada.append(magnitude_spectrum.copy())
            frames_transformada.append(power_spectrum.copy())
            kx_freqs.append(kx)
            kz_freqs.append(kz)

        # Avanzar en el tiempo
        psi_nm1 = psi_n.copy()
        psi_n = psi_np1.copy()

    return frames, E_history, frames_transformada, kx_freqs[0], kz_freqs[0], amplitud_en_coordenada

def calcular_energia_total(psi_n, psi_nm1, dx, dz, dt, g):
    """
    Calcula la energía total del sistema 
    
    Para ecuaciones de aguas poco profundas:
    E = E_cinetica + E_potencial
    E_cinetica = (1/2) * rho * H * |u|^2
    E_potencial = (1/2) * rho * g * eta^2
    
    donde u es la velocidad y eta es la elevación de la superficie
    """
    global H
    global zoom
    global X, Z, x

    # Crear copias locales para no modificar los argumentos originales
    psi_n_calc = psi_n
    psi_nm1_calc = psi_nm1
    H_calc = H
    
    if zoom:
        xmin, xmax = -10, 10
        
        # Verificar que x existe y tiene elementos
        if x is not None and len(x) > 0:
            zona_x = (x >= xmin) & (x <= xmax)
            zona_idx = np.where(zona_x)[0]
            
            # Verificar que hay índices válidos
            if len(zona_idx) > 0:
                # Usar slicing seguro para psi y H
                psi_n_calc = psi_n[zona_idx, :]
                psi_nm1_calc = psi_nm1[zona_idx, :]
                H_calc = H[zona_idx, :]  # ¡CRÍTICO! También recortar H
            else:
                print(f"Advertencia: No hay puntos en el rango [{xmin}, {xmax}]")
                print(f"Rango de x: [{x.min()}, {x.max()}]")
        else:
            print("Error: La variable global 'x' no está definida correctamente")
        
    # 1. Derivada temporal (velocidad vertical relacionada)
    Phi_t = (psi_n_calc - psi_nm1_calc) / dt
    
    # 2. Gradientes espaciales con manejo de bordes
    Phi_x = np.gradient(psi_n_calc, dx, axis=0)
    Phi_z = np.gradient(psi_n_calc, dz, axis=1)
    
    # 3. Energía cinética: (1/2) * H * (∂Φ/∂t)^2
    E_cinetica = 0.5 * H_calc * Phi_t**2
    
    # 4. Energía potencial: (1/2) * g * H * |∇Φ|^2
    E_potencial = 0.5 * g * H_calc * (Phi_x**2 + Phi_z**2)
    
    # 5. Energía total integrada sobre el dominio
    energy_density = E_cinetica + E_potencial
    E_total = np.sum(energy_density) * dx * dz
    
    return E_total, E_cinetica, E_potencial, energy_density