import datetime
import os
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.animation import PillowWriter
import numpy as np
from Source.calculos import calcular_canal

frames = None
E_history = None
dt = None
X = None
Z = None
frames_transformada = None
kx_freqs = None
kz_freqs = None
zoom = False
amplitud_en_coordenada = None

root = tk.Tk()
root.title("Interfaz")
root.geometry("1400x700")
root.iconbitmap(r"trappedmodes\Source\Logo3.ico")
root.config(bg='#bababa')

# Variables globales
Dimension = tk.StringVar(value="2D") #Por defecto 2D

# Variables
Longitud = tk.StringVar(value="10.0")
Ancho = tk.StringVar(value="2.0")
Profundidad = tk.StringVar(value="1.0")
Tiempo = tk.StringVar(value="10.0")
tipo_entrada = tk.StringVar(value="Water Drop")
tipo_perturbacion = tk.StringVar(value="Fondo Plano")
frontera_lat_izq = tk.StringVar(value="Infinite")
frontera_lat_der = tk.StringVar(value="Infinite")
frontera_lat_sup = tk.StringVar(value="Neumann")
frontera_lat_inf = tk.StringVar(value="Neumann")
Coordenada_x = tk.StringVar(value="0.0")
Coordenada_y = tk.StringVar(value="0.0")
Tipo_control = tk.StringVar(value="Deshabilitado")
Ganancia = tk.StringVar(value="0.0")
Frecuencia = tk.StringVar(value="")
Dato1 = tk.StringVar(value="")
Dato2 = tk.StringVar(value="")
Dato3 = tk.StringVar(value="")

# Variables para sliders
res_x = tk.IntVar(value=301)
res_z = tk.IntVar(value=201)

# Frame principal
main_frame = tk.Frame(root, bg='lightgray', relief=tk.RIDGE, borderwidth=2)
main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Título
title = tk.Label(main_frame, text="Wave Simulator for Trapped Mode Detection", font=('Arial', 13, 'bold'), bg='lightgray')
title.pack(pady=0.01)

# Frame contenedor principal
content_frame = tk.Frame(main_frame, bg='white', relief=tk.RIDGE, borderwidth=2)
content_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Configurar grid del content_frame
content_frame.columnconfigure(0, weight=0, minsize=280)  # Panel izquierdo
content_frame.columnconfigure(1, weight=1)  # Panel central (gráficas)
content_frame.columnconfigure(2, weight=0, minsize=280)  # Panel derecho
content_frame.rowconfigure(0, weight=1)

# ========================== PANEL IZQUIERDO ==========================
left_panel = tk.Frame(content_frame, bg='white', relief=tk.RIDGE, borderwidth=2)
left_panel.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)

# Frame para las 4 primeras entradas alineadas
inputs_frame = tk.Frame(left_panel, bg='white')
inputs_frame.pack(fill=tk.X, padx=5, pady=5)

# Configurar grid para alineación
inputs_frame.columnconfigure(0, weight=1, uniform="cols")
inputs_frame.columnconfigure(1, weight=1, uniform="cols")

tk.Label(left_panel, text="Simulation Parameters", bg='white', font=('Arial', 9, 'bold')).pack(pady=(0,5))

# Longitud del Canal
tk.Label(inputs_frame, text="Channel Length [m]", bg='white', font=('Arial', 8)).grid(row=0, column=0, padx=2, pady=2)
tk.Entry(inputs_frame, textvariable=Longitud, width=12, relief=tk.RIDGE, borderwidth=2, justify='center').grid(row=1, column=0, padx=2, pady=2)

# Ancho del Canal
tk.Label(inputs_frame, text="Channel Width [m]", bg='white', font=('Arial', 8)).grid(row=0, column=1, padx=2, pady=2)
tk.Entry(inputs_frame, textvariable=Ancho, width=12, relief=tk.RIDGE, borderwidth=2, justify='center').grid(row=1, column=1, padx=2, pady=2)

# Profundidad H
tk.Label(inputs_frame, text="Depth h [m]", bg='white', font=('Arial', 8)).grid(row=2, column=0, padx=2, pady=2)
tk.Entry(inputs_frame, textvariable=Profundidad, width=12, relief=tk.RIDGE, borderwidth=2, justify='center').grid(row=3, column=0, padx=2, pady=2)

# Tiempo de Simulación
tk.Label(inputs_frame, text="Simulation Time [s]", bg='white', font=('Arial', 8)).grid(row=2, column=1, padx=2, pady=2)
tk.Entry(inputs_frame, textvariable=Tiempo, width=12, relief=tk.RIDGE, borderwidth=2, justify='center').grid(row=3, column=1, padx=2, pady=2)


# Tipo de entrada y Tipo de Perturbación (misma fila)
row3_frame = tk.Frame(left_panel, bg='white')
row3_frame.pack(fill=tk.X, padx=5, pady=5)

col1_frame = tk.Frame(row3_frame, bg='white')
col1_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col1_frame, text="Input Type", bg='white', font=('Arial', 8)).pack()
combo1 = ttk.Combobox(col1_frame, textvariable=tipo_entrada, 
                      values=["Water Drop", "Gaussian", "Sinusoidal", "None"], 
                      state='readonly', width=13, font=('Arial', 8))
combo1.pack()

col2_frame = tk.Frame(row3_frame, bg='white')
col2_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col2_frame, text="Perturbation Type", bg='white', font=('Arial', 8)).pack()
combo2 = ttk.Combobox(col2_frame, textvariable=tipo_perturbacion, 
                      values=["Flat Bottom", "Rectangular", "Smooth Crest", "Gaussian"], 
                      state='readonly', width=12, font=('Arial', 8))
combo2.pack()

# Condiciones de frontera
tk.Label(left_panel, text="Boundary Conditions", bg='white', font=('Arial', 9, 'bold')).pack(pady=(10,5))

# Lateral Izq. y Lateral Der.
row4_frame = tk.Frame(left_panel, bg='white')
row4_frame.pack(fill=tk.X, padx=5, pady=2)

valores_frontera = ["Neumann", "Dirichlet", "Open", "Absorbing", "Absorbing2"]
valores_frontera_x = ["Neumann", "Dirichlet", "Open", "Absorbing", "Absorbing2", "Infinite"]

col1_frame = tk.Frame(row4_frame, bg='white')
col1_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col1_frame, text="Left", bg='white', font=('Arial', 8)).pack()
combo3 = ttk.Combobox(col1_frame, textvariable=frontera_lat_izq, 
                      values=valores_frontera_x, 
                      state='readonly', width=12, font=('Arial', 8))
combo3.pack()

col2_frame = tk.Frame(row4_frame, bg='white')
col2_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col2_frame, text="Right", bg='white', font=('Arial', 8)).pack()
combo4 = ttk.Combobox(col2_frame, textvariable=frontera_lat_der, 
                      values=valores_frontera_x, 
                      state='readonly', width=12, font=('Arial', 8))
combo4.pack()

# Lateral Sup. y Lateral Inf.
row5_frame = tk.Frame(left_panel, bg='white')
row5_frame.pack(fill=tk.X, padx=5, pady=2)

col1_frame = tk.Frame(row5_frame, bg='white')
col1_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col1_frame, text="Top", bg='white', font=('Arial', 8)).pack()
combo5 = ttk.Combobox(col1_frame, textvariable=frontera_lat_sup, 
                      values=valores_frontera, 
                      state='readonly', width=12, font=('Arial', 8))
combo5.pack()

col2_frame = tk.Frame(row5_frame, bg='white')
col2_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col2_frame, text="Bottom", bg='white', font=('Arial', 8)).pack()
combo6 = ttk.Combobox(col2_frame, textvariable=frontera_lat_inf, 
                      values=valores_frontera, 
                      state='readonly', width=12, font=('Arial', 8))
combo6.pack()

tk.Label(left_panel, text="Measurement Coordinates", bg='white', font=('Arial', 9, 'bold')).pack(pady=(10,5))

# Coordenada x y Coordenada z
row6_frame = tk.Frame(left_panel, bg='white')
row6_frame.pack(fill=tk.X, padx=5, pady=5)

col1_frame = tk.Frame(row6_frame, bg='white')
col1_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col1_frame, text="x Coordinate", bg='white', font=('Arial', 8)).pack()
tk.Entry(col1_frame, textvariable=Coordenada_x, width=12, relief=tk.RIDGE, borderwidth=2, justify='center').pack()

col2_frame = tk.Frame(row6_frame, bg='white')
col2_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col2_frame, text="z Coordinate", bg='white', font=('Arial', 8)).pack()
tk.Entry(col2_frame, textvariable=Coordenada_y, width=12, relief=tk.RIDGE, borderwidth=2, justify='center').pack()

#Parametros de Control
tk.Label(left_panel, text="Control Parameters", bg='white', font=('Arial', 9, 'bold')).pack(pady=(10,5))

# Tipo de Control y Ganancia
row7_frame = tk.Frame(left_panel, bg='white')
row7_frame.pack(fill=tk.X, padx=5, pady=5)

col1_frame = tk.Frame(row7_frame, bg='white')
col1_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col1_frame, text="Control Type", bg='white', font=('Arial', 8)).pack()
combo7 = ttk.Combobox(col1_frame, textvariable=Tipo_control, 
                      values=["Disabled", "Velocity Feedback"], 
                      state='readonly', width=12, font=('Arial', 8))
combo7.pack()

col2_frame = tk.Frame(row7_frame, bg='white')
col2_frame.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
tk.Label(col2_frame, text="Gain k", bg='white', font=('Arial', 8)).pack()
tk.Entry(col2_frame, textvariable=Ganancia, width=12, relief=tk.RIDGE, borderwidth=2, justify='center').pack()


# ********************************** FUNCIONES **********************************
def export_3d(frames, dt, frames_transformada, kx_freqs, kz_freqs):
    global zoom, X, Z, fig1, fig2, fig3
    
    # ===== CONFIGURACIÓN INICIAL =====
    output_path_base = "trappedmodes/Output"
    date_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_folder = os.path.join(output_path_base, f"simulation_{date_str}")
    os.makedirs(output_folder, exist_ok=True)

    # Guardar parámetros de simulación
    params_path = os.path.join(output_folder, "simulation_parameters.txt")
    with open(params_path, 'w') as f:
        f.write("Simulation Parameters:\n")
        f.write(f"Channel Length [m]: {Longitud.get()}\n")
        f.write(f"Channel Width [m]: {Ancho.get()}\n")
        f.write(f"Depth h [m]: {Profundidad.get()}\n")
        f.write(f"Simulation Time [s]: {Tiempo.get()}\n")
        f.write(f"Input Type: {tipo_entrada.get()}\n")
        f.write(f"Perturbation Type: {tipo_perturbacion.get()}\n")
        f.write(f"Left Boundary: {frontera_lat_izq.get()}\n")
        f.write(f"Right Boundary: {frontera_lat_der.get()}\n")
        f.write(f"Top Boundary: {frontera_lat_sup.get()}\n")
        f.write(f"Bottom Boundary: {frontera_lat_inf.get()}\n")
        f.write(f"x Coordinate: {Coordenada_x.get()}\n")
        f.write(f"z Coordinate: {Coordenada_y.get()}\n")
        f.write(f"Control Type: {Tipo_control.get()}\n")
        f.write(f"Gain k: {Ganancia.get()}\n")
    
    print(f"Parameters saved to {params_path}")

    try:
        vmin, vmax = np.min(frames), np.max(frames)
        writer = PillowWriter(fps=10)
        
        # Preparar zoom si es necesario
        if zoom:
            xmin, xmax = -10, 10
            x = np.linspace(-float(Longitud.get()), float(Longitud.get()), int(res_x.get()))
            mask_x = (x >= xmin) & (x <= xmax)
            X_zoom, Z_zoom = X[mask_x, :], Z[mask_x, :]
        else:
            mask_x = slice(None)
            X_zoom, Z_zoom = X, Z

        # ===== G1: ANIMACIÓN 3D =====
        print("Exporting 3D animation...")
        fig1.clear()
        ax1 = fig1.add_subplot(111, projection='3d')

        def configure_axis_3d():
            ax1.set_xlabel('X (m)', fontsize=9, labelpad=8)
            ax1.set_ylabel('Z (m)', fontsize=9, labelpad=8)
            ax1.set_zlabel('ψ(x,z)', fontsize=9, labelpad=8)
            ax1.set_title("3D Wave Propagation", fontsize=11, pad=15, fontweight='bold')
            ax1.set_facecolor('#f0f8ff')
            ax1.grid(True, alpha=0.2)
            ax1.xaxis.pane.fill = False
            ax1.yaxis.pane.fill = False
            ax1.zaxis.pane.fill = False
            ax1.set_zlim(vmin, vmax)

        surf1 = ax1.plot_surface(
            X_zoom, Z_zoom, frames[0][mask_x, :],
            cmap='ocean', alpha=0.85, vmin=vmin, vmax=vmax,
            antialiased=False, edgecolor='none',
            rcount=25, ccount=25, shade=False
        )
        configure_axis_3d()
        fig1.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10, pad=0.1)

        def update_3d(i):
            ax1.clear()
            surf = ax1.plot_surface(
                X_zoom, Z_zoom, frames[i][mask_x, :],
                cmap='ocean', alpha=0.85, vmin=vmin, vmax=vmax,
                antialiased=False, edgecolor='none',
                rcount=25, ccount=25, shade=False
            )
            configure_axis_3d()
            ax1.set_title(f"3D Wave Propagation - Frame {i+1}/{len(frames)} ", 
                         fontsize=10, pad=15, fontweight='bold')
            return surf,

        ani = animation.FuncAnimation(fig1, update_3d, frames=len(frames),
                                     interval=150, blit=False, repeat=False)
        
        output_path_3d = os.path.join(output_folder, "3D_wave_propagation.gif")
        ani.save(output_path_3d, writer=writer)
        print(f"✓ 3D Animation saved to {output_path_3d}")

        # ===== G2: HISTORIAL DE ENERGÍA =====
        print("Exporting energy history...")
        E0 = E_history[0]
        #Si E0 es 0, evitar división por cero
        E_norm = np.array(E_history) / E0 if E0 != 0 else np.array(E_history)
        time = np.arange(len(E_history)) * dt
        
        fig2.clear()
        ax2 = fig2.add_subplot(111)
        ax2.plot(time, E_norm, color='orange', linewidth=1, label='Total Energy')
        ax2.set_title('Normalized Energy vs Time', fontsize=10)
        ax2.set_xlabel('Time Step [s]')
        ax2.set_ylabel('E/E₀ (Normalized Energy)')
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='upper right')
        fig2.subplots_adjust(left=0.15, bottom=0.18)
        
        energy_path = os.path.join(output_folder, "energy_history.png")
        fig2.savefig(energy_path, dpi=150)
        print(f"✓ Energy history saved to {energy_path}")

        # ===== G3: TRANSFORMADA DE FOURIER =====
        print("Exporting power spectrum...")
        fig3.clear()
        ax3 = fig3.add_subplot(111)
        
        power_spectrum = np.abs(frames_transformada[0])**2
        log_power = np.log10(power_spectrum.T + 1e-10)
        kx, kz = kx_freqs, kz_freqs
        k_max_x, k_max_z = kx.max() * 0.5, kz.max() * 0.5
        
        im3 = ax3.imshow(
            log_power, origin='lower', cmap='hot',
            extent=[kx.min(), kx.max(), kz.min(), kz.max()],
            aspect='auto', interpolation='bilinear',
            vmin=np.percentile(log_power, 5),
            vmax=np.percentile(log_power, 99.5)
        )
        
        ax3.axhline(y=0, color='white', linestyle='--', alpha=0.3, linewidth=0.8)
        ax3.axvline(x=0, color='white', linestyle='--', alpha=0.3, linewidth=0.8)
        ax3.set_xlim([-k_max_x, k_max_x])
        ax3.set_ylim([-k_max_z, k_max_z])
        ax3.set_title('Power Spectrum (2D FFT)', fontsize=12, fontweight='bold', pad=10)
        ax3.set_xlabel('$k_x$ (cycles/m)', fontsize=10)
        ax3.set_ylabel('$k_z$ (cycles/m)', fontsize=10)
        ax3.grid(True, alpha=0.2, linestyle=':', linewidth=0.5, color='white')
        fig3.subplots_adjust(left=0.15, bottom=0.18)

        def update_fft(frame_idx):
            power_spectrum = np.abs(frames_transformada[frame_idx])**2
            log_power = np.log10(power_spectrum.T + 1e-10)
            im3.set_array(log_power)
            ax3.set_title(f'Power Spectrum - Frame {frame_idx+1}/{len(frames_transformada)}',
                         fontsize=12, fontweight='bold', pad=10)
            return [im3]

        ani_fft = animation.FuncAnimation(fig3, update_fft, 
                                         frames=len(frames_transformada),
                                         interval=150, blit=False, repeat=False)
        
        fft_path = os.path.join(output_folder, "power_spectrum.gif")
        ani_fft.save(fft_path, writer=writer)
        print(f"✓ Power spectrum saved to {fft_path}")
        
        print(f"\n{'='*50}")
        print(f"Export completed successfully!")
        print(f"All files saved to: {output_folder}")
        print(f"{'='*50}\n")

    except Exception as e:
        print(f"❌ Error exporting 3D: {e}")
        import traceback
        traceback.print_exc()

def export_2d(frames, dt, frames_transformada, kx_freqs, kz_freqs):
    global zoom
    global fig1, fig2, fig3, fig4
    output_path_base = "trappedmodes/Output"
    # Crear folder con fecha y dentro la animación
    date_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_folder = os.path.join(output_path_base, f"simulation_{date_str}")
    os.makedirs(output_folder, exist_ok=True)

    # Salvar parametros en TXT
    params_path = os.path.join(output_folder, "simulation_parameters.txt")

    with open(params_path, 'w') as f:
        f.write("Simulation Parameters:\n")
        f.write(f"Channel Length [m]: {Longitud.get()}\n")
        f.write(f"Channel Width [m]: {Ancho.get()}\n")
        f.write(f"Depth h [m]: {Profundidad.get()}\n")
        f.write(f"Simulation Time [s]: {Tiempo.get()}\n")
        f.write(f"Input Type: {tipo_entrada.get()}\n")
        f.write(f"Perturbation Type: {tipo_perturbacion.get()}\n")
        f.write(f"Left Boundary: {frontera_lat_izq.get()}\n")
        f.write(f"Right Boundary: {frontera_lat_der.get()}\n")
        f.write(f"Top Boundary: {frontera_lat_sup.get()}\n")
        f.write(f"Bottom Boundary: {frontera_lat_inf.get()}\n")
        f.write(f"x Coordinate: {Coordenada_x.get()}\n")
        f.write(f"z Coordinate: {Coordenada_y.get()}\n")
        f.write(f"Control Type: {Tipo_control.get()}\n")
        f.write(f"Gain k: {Ganancia.get()}\n")

    try:
        vmin = np.min(frames)
        vmax = np.max(frames)
        
        # G1: Wave propagation
        fig1.clear()    
        ax1 = fig1.add_subplot(111)

        if zoom:
            xmin, xmax = -10, 10
            x = np.linspace( -float(Longitud.get()), float(Longitud.get()), int(res_x.get()))
            mask_x = (x >= xmin) & (x <= xmax)
            X_zoom = X[mask_x, :]
            Z_zoom = Z[mask_x, :]
            extent = [X_zoom.min(), X_zoom.max(), Z_zoom.min(), Z_zoom.max()]
        else:
            mask_x = slice(None)  # No filtrar
            extent = [-float(Longitud.get()), float(Longitud.get()), -float(Ancho.get()), float(Ancho.get())]

        # Imagen inicial
        im = ax1.imshow(
            frames[0][mask_x, :].T,
            origin='lower',
            cmap='ocean',
            extent=extent,
            aspect='auto',
            interpolation='bilinear',
            vmin=vmin, vmax=vmax
        )

        #ax1.set_title("Solution with constant depth h", fontsize=12, fontweight='bold', pad=10)
        
        ax1.set_title("Wave Propagation", fontsize=12, fontweight='bold', pad=10)
        ax1.set_xlabel('X (m)', fontsize=10)
        ax1.set_ylabel('Z (m)', fontsize=10)
        ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

        fig1.colorbar(im, ax=ax1, shrink=0.5, aspect=10)
        
        fig1.subplots_adjust(left=0.15, bottom=0.18)

        def update(frame_idx):
            im.set_array(frames[frame_idx][mask_x, :].T)
            #Titulo Frame actual
            ax1.set_title(f"Wave Propagation - Frame {frame_idx+1}/{len(frames)}", fontsize=12, fontweight='bold', pad=10)
            return [im]
        
        ani = animation.FuncAnimation(
            fig1,
            update,
            frames=len(frames),
            interval=150,   # ms
            blit=False,
            repeat=False
        )
        # Save GIF
        output_path = os.path.join(output_folder, "wave_propagation.gif")
        writer = PillowWriter(fps=10)
        ani.save(output_path, writer=writer)
        print(f"Animation saved to {output_path}")

        # Salvar G2: Energy history
        E0 = E_history[0]
        #Si E0 es 0, evitar división por cero
        E_norm = np.array(E_history) / E0 if E0 != 0 else np.array(E_history)
        time = np.arange(len(E_history)) * dt
        fig2.clear()
        ax2 = fig2.add_subplot(111)
        ax2.plot(time, E_norm, color='orange', linewidth=1, label='Total Energy')
        ax2.set_title('Normalized Energy vs Time', fontsize=10)
        ax2.set_xlabel('Time Step [s]')
        ax2.set_ylabel('E/E₀ (Normalized Energy)')
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='upper right')
        fig2.subplots_adjust(left=0.15, bottom=0.18)
        energy_path = os.path.join(output_folder, "energy_history.png")
        fig2.savefig(energy_path)
        print(f"Energy history saved to {energy_path}")

        # Salvar G3: Transformada de Fourier como GIF
        fig3.clear()
        ax3 = fig3.add_subplot(111)
        power_spectrum = np.abs(frames_transformada[0])**2
        log_power = np.log10(power_spectrum.T + 1e-10)
        kx = kx_freqs
        kz = kz_freqs
        k_max_x = kx.max() * 0.5
        k_max_z = kz.max() * 0.5
        im3 = ax3.imshow(
            log_power,
            origin='lower',
            cmap='hot',
            extent=[kx.min(), kx.max(), kz.min(), kz.max()],
            aspect='auto',
            interpolation='bilinear',
            vmin=np.percentile(log_power, 5),
            vmax=np.percentile(log_power, 99.5)
        )
        ax3.axhline(y=0, color='white', linestyle='--', alpha=0.3, linewidth=0.8)
        ax3.axvline(x=0, color='white', linestyle='--', alpha=0.3, linewidth=0.8)
        ax3.set_xlim([-k_max_x, k_max_x])
        ax3.set_ylim([-k_max_z, k_max_z])
        ax3.set_title('(2D FFT)', fontsize=10, fontweight='bold', pad=10)
        ax3.set_xlabel('$k_x$ (cycles/m)', fontsize=10)
        ax3.set_ylabel('$k_z$ (cycles/m)', fontsize=10)
        ax3.grid(True, alpha=0.2, linestyle=':', linewidth=0.5, color='white')
        fig3.colorbar(im3, ax=ax3, shrink=0.5, aspect=10)
        fig3.subplots_adjust(left=0.15, bottom=0.18)
        
        def update_fft(frame_idx):
            power_spectrum = np.abs(frames_transformada[frame_idx])**2
            log_power = np.log10(power_spectrum.T + 1e-10)
            im3.set_array(log_power)
            ax3.set_title(f'(2D FFT) - Frame {frame_idx+1}/{len(frames_transformada)}', 
                         fontsize=10, fontweight='bold', pad=10)
            return [im3]
        ani_fft = animation.FuncAnimation(
            fig3,
            update_fft,
            frames=len(frames_transformada),
            interval=150,
            blit=False,
            repeat=False
        )
        # Save GIF
        fft_path = os.path.join(output_folder, "power_spectrum.gif")
        ani_fft.save(fft_path, writer=writer)
        print(f"Power spectrum animation saved to {fft_path}")
        

    except Exception as e:
        print(f"Error al graficar en 2D: {e}")

def plot_2d(frames, dt, frames_transformada, kx_freqs, kz_freqs):
    global ax1, ax2, ax3, ax4
    global im, ani, im3, ani_fft
    global E_history
    global zoom

    # Detener animaciones previas de forma segura
    try:
        if 'ani' in globals() and ani is not None and hasattr(ani, 'event_source'):
            ani.event_source.stop()
    except:
        pass
    try:
        if 'ani_fft' in globals() and ani_fft is not None and hasattr(ani_fft, 'event_source'):
            ani_fft.event_source.stop()
    except:
        pass
    

    try:
        vmin = np.min(frames)
        vmax = np.max(frames)
        
        # G1: Wave propagation
        fig1.clear()    
        ax1 = fig1.add_subplot(111)

        if zoom:
            xmin, xmax = -10, 10
            x = np.linspace( -float(Longitud.get()), float(Longitud.get()), int(res_x.get()))
            mask_x = (x >= xmin) & (x <= xmax)
            X_zoom = X[mask_x, :]
            Z_zoom = Z[mask_x, :]
            extent = [X_zoom.min(), X_zoom.max(), Z_zoom.min(), Z_zoom.max()]
        else:
            mask_x = slice(None)  # No filtrar
            extent = [-float(Longitud.get()), float(Longitud.get()), -float(Ancho.get()), float(Ancho.get())]

        # Imagen inicial
        im = ax1.imshow(
            frames[0][mask_x, :].T,
            origin='lower',
            cmap='ocean',
            extent=extent,
            aspect='auto',
            interpolation='bilinear',
            vmin=vmin, vmax=vmax
        )

        #ax1.set_title("Solution with constant depth h", fontsize=12, fontweight='bold', pad=10)
        
        ax1.set_title("Wave Propagation", fontsize=12, fontweight='bold', pad=10)
        ax1.set_xlabel('X (m)', fontsize=10)
        ax1.set_ylabel('Z (m)', fontsize=10)
        ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
        #Graficar coordenadas de medición
        x_med = float(Coordenada_x.get())
        z_med = float(Coordenada_y.get())
        ax1.plot(x_med, z_med, marker='o', color='red', markersize=1, label='Measurement Point')

        fig1.colorbar(im, ax=ax1, shrink=0.5, aspect=10)
        #cbar1.set_label('Surface Height (m)', rotation=270, labelpad=20, fontsize=9)

        fig1.subplots_adjust(left=0.15, bottom=0.18)

        def update(frame_idx):
            im.set_array(frames[frame_idx][mask_x, :].T)
            #Titulo Frame actual con tiempo = frame_idx * (total_time / total_frames)
            ax1.set_title(f"Frame {frame_idx+1}/{len(frames)}", fontsize=10, fontweight='bold', pad=10)
            return [im]
        
        ani = animation.FuncAnimation(
            fig1,
            update,
            frames=len(frames),
            interval=150,   # ms
            blit=False,
            repeat=False
        )

        canvas1.draw()
            
        # G2: Energy history
        # Graficar energía normalizada vs tiempo
        E0 = E_history[0]
        #Si E0 es 0, evitar división por cero
        E_norm = np.array(E_history) / E0 if E0 != 0 else np.array(E_history)
        time = np.arange(len(E_history)) * dt
        fig2.clear()
        ax2 = fig2.add_subplot(111)
        ax2.plot(time, E_norm, color='orange', linewidth=1, label='Total Energy')
        ax2.set_title('Normalized Energy vs Time', fontsize=10)
        ax2.set_xlabel('Time Step [s]')
        ax2.set_ylabel('E/E₀ (Normalized Energy)')
        ax2.grid(True, alpha=0.3) 
        ax2.legend(loc='upper right')
        fig2.subplots_adjust(left=0.15, bottom=0.18)
        canvas2.draw()
            
        # G3: Transformada de Fourier 
        fig3.clear()
        ax3 = fig3.add_subplot(111)
        power_spectrum = np.abs(frames_transformada[0])**2
        
        # Aplicar log y normalizar para mejor contraste
        log_power = np.log10(power_spectrum.T + 1e-10)
        
        # Límites de frecuencia (usar las frecuencias reales)
        kx = kx_freqs
        kz = kz_freqs
        
        # Encontrar límites razonables (excluyendo extremos con poco contenido)
        k_max_x = kx.max() * 0.5  # Mostrar hasta la mitad de la frecuencia máxima
        k_max_z = kz.max() * 0.5
        
        im3 = ax3.imshow(
            log_power,
            origin='lower',
            cmap='hot',  # 'hot', 'jet', o 'plasma' son buenos para espectros
            extent=[kx.min(), kx.max(), kz.min(), kz.max()],
            aspect='auto',
            interpolation='bilinear',
            vmin=np.percentile(log_power, 5),   # Mejorar contraste
            vmax=np.percentile(log_power, 99.5)
        )
        
        # Añadir líneas de referencia en k=0
        ax3.axhline(y=0, color='white', linestyle='--', alpha=0.3, linewidth=0.8)
        ax3.axvline(x=0, color='white', linestyle='--', alpha=0.3, linewidth=0.8)
        
        # Zoom a la región de interés (frecuencias bajas)
        ax3.set_xlim([-k_max_x, k_max_x])
        ax3.set_ylim([-k_max_z, k_max_z])
        
        ax3.set_title('(2D FFT)', fontsize=12, fontweight='bold', pad=10)
        ax3.set_xlabel('$k_x$ (cycles/m)', fontsize=10)
        ax3.set_ylabel('$k_z$ (cycles/m)', fontsize=10)
        ax3.grid(True, alpha=0.2, linestyle=':', linewidth=0.5, color='white')
        
        fig3.colorbar(im3, ax=ax3, shrink=0.5, aspect=10)
        #cbar3.set_label('$\\log_{10}(|\\hat{\\psi}|^2)$', rotation=270, labelpad=20, fontsize=9)
        
        fig3.subplots_adjust(left=0.15, bottom=0.18)

        def update_fft(frame_idx):
            # Actualizar con espectro de potencia
            power_spectrum = np.abs(frames_transformada[frame_idx])**2
            log_power = np.log10(power_spectrum.T + 1e-10)
            
            im3.set_array(log_power)
            
            # Actualizar límites de color para mantener contraste
            #im3.set_clim(vmin=np.percentile(log_power, 5), 
            #            vmax=np.percentile(log_power, 99.5))
            
            ax3.set_title(f'(2D FFT) - Frame {frame_idx+1}/{len(frames_transformada)}', 
                         fontsize=10, fontweight='bold', pad=10)
            return [im3]
        
        ani_fft = animation.FuncAnimation(
            fig3,
            update_fft,
            frames=len(frames_transformada),
            interval=150,
            blit=False,
            repeat=False
        )

        canvas3.draw()
            
        # G4: Amplitud en punto específico (x,z)
        fig4.clear()
        ax4 = fig4.add_subplot(111)
        global amplitud_en_coordenada
        ax4.plot(time, amplitud_en_coordenada, color='green', linewidth=1)
        ax4.set_title(f'Amplitude at (x={Coordenada_x.get()} m, z={Coordenada_y.get()} m) vs Time', fontsize=10)
        ax4.set_xlabel('Time Step [s]')
        ax4.set_ylabel('Amplitude ψ(x,z)')
        ax4.grid(True, alpha=0.3)   
        fig4.subplots_adjust(left=0.15, bottom=0.18)
        canvas4.draw()

    except Exception as e:
        print(f"Error al graficar en 2D: {e}")
    
def plot_3d(frames, X, Z, dt, frames_transformada, kx_freqs, kz_freqs):
    """Grafica funciones en 3D"""
    # Crear malla para gráficas 3D
    vmin = np.min(frames)
    vmax = np.max(frames)

    global ax1, ax2, surf1, ani
    global im3, ani_fft
    global E_history
    global zoom

    # Detener animaciones previas de forma segura
    try:
        if 'ani' in globals() and ani is not None and hasattr(ani, 'event_source'):
            ani.event_source.stop()
    except:
        pass
    try:
        if 'ani_fft' in globals() and ani_fft is not None and hasattr(ani_fft, 'event_source'):
            ani_fft.event_source.stop()
    except:
        pass

    fig1.clear()
    ax1 = fig1.add_subplot(111, projection='3d')  

    if zoom:
            xmin, xmax = -10, 10
            x = np.linspace( -float(Longitud.get()), float(Longitud.get()), int(res_x.get()))
            mask_x = (x >= xmin) & (x <= xmax)
            X_zoom = X[mask_x, :]
            Z_zoom = Z[mask_x, :]
    else:
        mask_x = slice(None)
        X_zoom = X
        Z_zoom = Z

    def configure_axis():
        ax1.set_xlabel('X (m)', fontsize=9, labelpad=8)
        ax1.set_ylabel('Z (m)', fontsize=9, labelpad=8)
        ax1.set_zlabel('ψ(x,z)', fontsize=9, labelpad=8)
        ax1.set_title("3D Wave Propagation", fontsize=11, pad=15, fontweight='bold')
        ax1.set_facecolor('#f0f8ff')
        ax1.grid(True, alpha=0.2)
        ax1.xaxis.pane.fill = False
        ax1.yaxis.pane.fill = False
        ax1.zaxis.pane.fill = False
        ax1.set_zlim(vmin, vmax)
        #ax1.set_box_aspect((2.5, 1.2, 0.4))

    # Surface inicial
    surf1 = ax1.plot_surface(
        X_zoom, Z_zoom, frames[0][mask_x, :],
        cmap='ocean',
        alpha=0.85,
        vmin=vmin, vmax=vmax,
        antialiased=False,
        edgecolor='none',
        rcount=25, ccount=25,
        shade=False
    )

    configure_axis()

    def update(i):
        ax1.clear()
        surf = ax1.plot_surface(
            X_zoom, Z_zoom, frames[i][mask_x, :],
            cmap='ocean',
            alpha=0.85,
            vmin=vmin, vmax=vmax,
            antialiased=False,
            edgecolor='none',
            rcount=25, ccount=25,
            shade=False
        )
        configure_axis()
        ax1.set_zlim(vmin, vmax)
        ax1.set_title(f"3D - Frame {i+1}/{len(frames)} ", fontsize=11, pad=15, fontweight='bold')
        return surf,


    ani = animation.FuncAnimation(
        fig1,
        update,
        frames=len(frames),
        interval=150, 
        blit=False,
        repeat=False
    )
    
    canvas1.draw()   
        
    # G2: Energy history
    # Graficar energía normalizada vs tiempo
    E0 = E_history[0]
    #Si E0 es 0, evitar división por cero
    E_norm = np.array(E_history) / E0 if E0 != 0 else np.array(E_history)
    time = np.arange(len(E_history)) * dt
    fig2.clear()
    ax2 = fig2.add_subplot(111)
    ax2.plot(time, E_norm, color='orange', linewidth=1, label='Total Energy')
    ax2.set_title('Normalized Energy vs Time', fontsize=10)
    ax2.set_xlabel('Time Step [s]')
    ax2.set_ylabel('E/E₀ (Normalized Energy)')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')
    fig2.subplots_adjust(left=0.15, bottom=0.18)
    canvas2.draw()
        
    # G3: Transformada de Fourier 
    fig3.clear()
    ax3 = fig3.add_subplot(111)
    power_spectrum = np.abs(frames_transformada[0])**2
    
    # Aplicar log y normalizar para mejor contraste
    log_power = np.log10(power_spectrum.T + 1e-10)
    
    # Límites de frecuencia (usar las frecuencias reales)
    kx = kx_freqs
    kz = kz_freqs
    
    # Encontrar límites razonables (excluyendo extremos con poco contenido)
    k_max_x = kx.max() * 0.5  # Mostrar hasta la mitad de la frecuencia máxima
    k_max_z = kz.max() * 0.5
    
    im3 = ax3.imshow(
        log_power,
        origin='lower',
        cmap='hot',  # 'hot', 'jet', o 'plasma' son buenos para espectros
        extent=[kx.min(), kx.max(), kz.min(), kz.max()],
        aspect='auto',
        interpolation='bilinear',
        vmin=np.percentile(log_power, 5),   # Mejorar contraste
        vmax=np.percentile(log_power, 99.5)
    )
    
    # Añadir líneas de referencia en k=0
    ax3.axhline(y=0, color='white', linestyle='--', alpha=0.3, linewidth=0.8)
    ax3.axvline(x=0, color='white', linestyle='--', alpha=0.3, linewidth=0.8)
    
    # Zoom a la región de interés (frecuencias bajas)
    ax3.set_xlim([-k_max_x, k_max_x])
    ax3.set_ylim([-k_max_z, k_max_z])
    
    ax3.set_title('Power Spectrum (2D FFT)', fontsize=10, fontweight='bold', pad=10)
    ax3.set_xlabel('$k_x$ (cycles/m)', fontsize=10)
    ax3.set_ylabel('$k_z$ (cycles/m)', fontsize=10)
    ax3.grid(True, alpha=0.2, linestyle=':', linewidth=0.5, color='white')
    
    fig3.colorbar(im3, ax=ax3, shrink=0.5, aspect=10)
    #cbar3.set_label('$\\log_{10}(|\\hat{\\psi}|^2)$', rotation=270, labelpad=20, fontsize=9)
    
    fig3.subplots_adjust(left=0.15, bottom=0.18)

    def update_fft(frame_idx):
        # Actualizar con espectro de potencia
        power_spectrum = np.abs(frames_transformada[frame_idx])**2
        log_power = np.log10(power_spectrum.T + 1e-10)
        
        im3.set_array(log_power)
        
        # Actualizar límites de color para mantener contraste
        #im3.set_clim(vmin=np.percentile(log_power, 5), 
        #            vmax=np.percentile(log_power, 99.5))
        
        ax3.set_title(f'Power Spectrum (2D FFT) - Frame {frame_idx+1}/{len(frames_transformada)}', 
                        fontsize=12, fontweight='bold', pad=10)
        return [im3]
    
    ani_fft = animation.FuncAnimation(
        fig3,
        update_fft,
        frames=len(frames_transformada),
        interval=150,
        blit=False,
        repeat=False
    )

    canvas3.draw()
        
    # G4: Amplitud en punto específico (x,z)
    fig4.clear()
    ax4 = fig4.add_subplot(111)
    global amplitud_en_coordenada
    ax4.plot(time, amplitud_en_coordenada, color='green', linewidth=1)
    ax4.set_title(f'Amplitude at (x={Coordenada_x.get()} m, z={Coordenada_y.get()} m) vs Time', fontsize=10)
    ax4.set_xlabel('Time Step [s]')
    ax4.set_ylabel('Amplitude ψ(x,z)')
    ax4.grid(True, alpha=0.3)   
    fig4.subplots_adjust(left=0.15, bottom=0.18)
    canvas4.draw()

def Calcular():      
    global frames
    global X, Z
    global E_history
    global dt
    global frames_transformada
    global kx_freqs, kz_freqs
    global amplitud_en_coordenada

    print("=== COMPUTATION ===")
    Dato3.set("Computing...")
    # Color
    Dato3_label = tk.Label(left_panel, textvariable=Dato3, bg='white', font=('Arial', 8, 'italic'), fg='blue')
    Dato3_label.pack(pady=5)
    root.update_idletasks()

    try:
        lon = float(Longitud.get())
        ancho = float(Ancho.get())
        res_x_val = res_x.get()
        res_z_val = res_z.get()
        h = float(Profundidad.get()) # if Profundidad.get() else 1.0
        tiempo = float(Tiempo.get()) # if Tiempo.get() else 5.0
        # condicion_inicial igual a tipo_entrada
        condicion_inicial = tipo_entrada.get()
        condicion_frontera = [frontera_lat_izq.get(), frontera_lat_der.get(), frontera_lat_sup.get(), frontera_lat_inf.get()]
        if condicion_frontera[0] == "Infinite" or condicion_frontera[1] == "Infinite":
            #print("Condición de frontera 'Infinite' seleccionada.")
            res_x_val = 601  # Ajustar resolución en x para frontera infinita
            lon = 100.0    # Ajustar longitud del canal para frontera infinita
            Longitud.set(str(lon))
            res_x.set(res_x_val)
            global zoom
            zoom = True
            frontera_lat_der.set("Infinite")
            frontera_lat_izq.set("Infinite")

        tipo_perturbacion_val = tipo_perturbacion.get()
        control_tipo = Tipo_control.get()
        ganancia_val = float(Ganancia.get()) if Ganancia.get() else 0
        coordenadas = (float(Coordenada_x.get()), float(Coordenada_y.get())) if Coordenada_x.get() and Coordenada_y.get() else (0.0, 0.0)

        frames, X, Z, E_history, dt, frames_transformada, kx_freqs, kz_freqs, amplitud_en_coordenada = calcular_canal(
            lon, ancho, Nx=res_x_val, Nz=res_z_val, h=h, 
            T=tiempo, condicion_inicial=condicion_inicial, 
            condicion_frontera=condicion_frontera, 
            tipo_perturbacion=tipo_perturbacion_val,
            control_tipo=control_tipo,
            ganancia_val=ganancia_val,coordenadas=coordenadas)
        print("Computation completed.")
        print("===============\n")
        Dato3.set("Computation completed.")
        # Color
        Dato3_label = tk.Label(left_panel, textvariable=Dato3, bg='white', font=('Arial', 8, 'italic'), fg='blue')
        Dato3_label.pack(pady=5)
    except ValueError:
        print("Error: Entrada inválida. Por favor, ingrese valores numéricos correctos.")
        Dato3.set("Error en entrada de datos.")
    
def Graficar():
    print("=== Plotting... ===")
    Dato3.set("Plotting...")
    Dato3_label = tk.Label(left_panel, textvariable=Dato3, bg='white', font=('Arial', 8, 'italic'), fg='blue')
    Dato3_label.pack(pady=5)
    root.update_idletasks()

    # Limpiar gráficas anteriores
    global fig1, fig2, fig3, fig4
    fig1.clf()
    fig2.clf()
    fig3.clf()
    fig4.clf()

    # Grafica las funciones en 2D o 3D
    if str(Dimension.get()) == "2D":
        plot_2d(frames, dt, frames_transformada, kx_freqs, kz_freqs)
    else:
        plot_3d(frames, X, Z, dt, frames_transformada, kx_freqs, kz_freqs)
    
    print("=== Plot completed. ===")
    print("===============\n")
    Dato3.set("Plot completed.")
    Dato3_label = tk.Label(left_panel, textvariable=Dato3, bg='white', font=('Arial', 8, 'italic'), fg='blue')
    Dato3_label.pack(pady=5)
    root.update_idletasks()

def Exportar():
    print("=== Exporting... ===")
    Dato3.set("Exporting...")
    Dato3_label = tk.Label(left_panel, textvariable=Dato3, bg='white', font=('Arial', 8, 'italic'), fg='blue')
    Dato3_label.pack(pady=5)
    root.update_idletasks()

    if str(Dimension.get()) == "2D":
        export_2d(frames, dt, frames_transformada, kx_freqs, kz_freqs)
    else:
        export_3d(frames, dt, frames_transformada, kx_freqs, kz_freqs)
    print("Export completed.")
    print("===============\n")
    Dato3.set("Export completed.")
    # Color
    Dato3_label = tk.Label(left_panel, textvariable=Dato3, bg='white', font=('Arial', 8, 'italic'), fg='blue')
    Dato3_label.pack(pady=5)
    root.update_idletasks()

def set_dimension(mode):
    # Cambia entre modo 2D y 3D
    Dimension.set(mode)
        
    if str(Dimension.get()) == "2D":
        btn_2d.config(relief=tk.SUNKEN, bg='#20AD4E', fg='white', borderwidth=3)
        btn_3d.config(relief=tk.RAISED, bg='lightgray', fg='black', borderwidth=2)
        print("Modo 2D activado")
    else:
        btn_3d.config(relief=tk.SUNKEN, bg='#20AD4E', fg='white', borderwidth=3)
        btn_2d.config(relief=tk.RAISED, bg='lightgray', fg='black', borderwidth=2)
        print("Modo 3D activado")

def stop_animations():
    global ani, ani_fft
    try:
        if 'ani' in globals() and ani is not None and hasattr(ani, 'event_source'):
            ani.event_source.stop()
            print("Animation 1 stopped.")
    except Exception as e:
        print(f"Error stopping animation 1: {e}")
    try:
        if 'ani_fft' in globals() and ani_fft is not None and hasattr(ani_fft, 'event_source'):
            ani_fft.event_source.stop()
            print("Animation 2 stopped.")
    except Exception as e:
        print(f"Error stopping animation 2: {e}")
# ********************************** FIN DE FUNCIONES **********************************


# ********************************** ESTRUCTURA GRÁFICAS **********************************

# ========================== PANEL CENTRAL (GRÁFICAS) ==========================
center_panel = tk.Frame(content_frame, bg='white')
center_panel.grid(row=0, column=1, sticky="nsew", padx=5, pady=5)

# Configurar grid 2x2
center_panel.rowconfigure(0, weight=1)
center_panel.rowconfigure(1, weight=1)
center_panel.columnconfigure(0, weight=1)
center_panel.columnconfigure(1, weight=1)

# G1
g1_frame = tk.Frame(center_panel, relief=tk.RIDGE, borderwidth=2)
g1_frame.grid(row=0, column=0, sticky="nsew", padx=3, pady=3)
tk.Label(g1_frame, text="Wave Propagation", bg='white', font=('Arial', 10)).pack()
fig1 = Figure(figsize=(4, 2.5), facecolor='white')
ax1 = fig1.add_subplot(111)
canvas1 = FigureCanvasTkAgg(fig1, g1_frame)
canvas1.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# G3
g3_frame = tk.Frame(center_panel, relief=tk.RIDGE, borderwidth=2)
g3_frame.grid(row=0, column=1, sticky="nsew", padx=3, pady=3)
tk.Label(g3_frame, text="Power Spectrum", bg='white', font=('Arial', 10)).pack()
fig3 = Figure(figsize=(4, 2.5), facecolor='white')
ax3 = fig3.add_subplot(111)
canvas3 = FigureCanvasTkAgg(fig3, g3_frame)
canvas3.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# G2
g2_frame = tk.Frame(center_panel, relief=tk.RIDGE, borderwidth=2)
g2_frame.grid(row=1, column=0, sticky="nsew", padx=3, pady=3)
tk.Label(g2_frame, text="Normalized Energy", bg='white', font=('Arial', 10)).pack()
fig2 = Figure(figsize=(4, 2.5), facecolor='white')
ax2 = fig2.add_subplot(111)
canvas2 = FigureCanvasTkAgg(fig2, g2_frame)
canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# G4
g4_frame = tk.Frame(center_panel, relief=tk.RIDGE, borderwidth=2)
g4_frame.grid(row=1, column=1, sticky="nsew", padx=3, pady=3)
tk.Label(g4_frame, text="Amplitude", bg='white', font=('Arial', 10)).pack()
fig4 = Figure(figsize=(4, 2.5), facecolor='white')
ax4 = fig4.add_subplot(111)
canvas4 = FigureCanvasTkAgg(fig4, g4_frame)
canvas4.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# ********************************** FIN DE ESTRUCTURA GRÁFICAS **********************************


# ========================== PANEL DERECHO ==========================
right_panel = tk.Frame(content_frame, bg='white', relief=tk.RIDGE, borderwidth=2)
right_panel.grid(row=0, column=2, sticky="nsew", padx=5, pady=5)

# Resolución en x con slider personalizado
tk.Label(right_panel, text="x Resolution", bg='white', font=('Arial', 9)).pack(pady=(10,5))
slider_frame_x = tk.Frame(right_panel, bg='white')
slider_frame_x.pack(padx=10)

slider_x = tk.Scale(slider_frame_x, from_=1, to=601, orient=tk.HORIZONTAL, 
                    variable=res_x, length=200, resolution=1,
                    bg='white', fg='black', troughcolor='#e0e0e0',
                    highlightthickness=0, bd=1, relief=tk.FLAT,
                    sliderlength=20, width=12, font=('Arial', 8))
slider_x.pack()

# Resolución en z con slider personalizado
tk.Label(right_panel, text="z Resolution", bg='white', font=('Arial', 9)).pack(pady=(10,5))
slider_frame_z = tk.Frame(right_panel, bg='white')
slider_frame_z.pack(padx=10)

slider_z = tk.Scale(slider_frame_z, from_=1, to=250, orient=tk.HORIZONTAL, 
                    variable=res_z, length=200, resolution=1,
                    bg='white', fg='black', troughcolor='#e0e0e0',
                    highlightthickness=0, bd=1, relief=tk.FLAT,
                    sliderlength=20, width=12, font=('Arial', 8))
slider_z.pack()

# ********************************** BOTONES **********************************

# Botones en dos filas
buttons_frame = tk.Frame(right_panel, bg='white')
buttons_frame.pack(pady=15)

# Primera fila de botones
btn_row1 = tk.Frame(buttons_frame, bg='white')
btn_row1.pack(pady=3)

btn_calcular = tk.Button(btn_row1, text="Compute", bg='#58DB83',
    relief=tk.RAISED, borderwidth=2, width=12,
    command= Calcular, font=('Arial', 9, 'bold'))
btn_calcular.pack(side=tk.LEFT, padx=3)

btn_stop = tk.Button(btn_row1, text="Stop", bg='#58DB83',
    relief=tk.RAISED, borderwidth=2, width=12,
    command=stop_animations, font=('Arial', 9, 'bold'))
btn_stop.pack(side=tk.LEFT, padx=3)

# Segunda fila de botones
btn_row2 = tk.Frame(buttons_frame, bg='white')
btn_row2.pack(pady=3)

btn_graficar = tk.Button(btn_row2, text="Plot", bg='#58DB83',
    relief=tk.RAISED, borderwidth=2, width=12,
    command= Graficar, font=('Arial', 9, 'bold'))
btn_graficar.pack(side=tk.LEFT, padx=3)

btn_exportar = tk.Button(btn_row2, text="Export", bg='#58DB83',
    relief=tk.RAISED, borderwidth=2, width=12,
    command=Exportar, font=('Arial', 9, 'bold'))
btn_exportar.pack(side=tk.LEFT, padx=3)

# Botones 2D y 3D
dimension_buttons = tk.Frame(right_panel, bg='white')
dimension_buttons.pack(pady=10)

btn_2d = tk.Button(dimension_buttons, text="2D", bg='white',
    relief=tk.RAISED, borderwidth=2, width=12,
    command=lambda: set_dimension("2D"), font=('Arial', 9, 'bold'))
btn_2d.pack(side=tk.LEFT, padx=3)

btn_3d = tk.Button(dimension_buttons, text="3D", bg='white',
    relief=tk.RAISED, borderwidth=2, width=12,
    command=lambda: set_dimension("3D"), font=('Arial', 9, 'bold'))
btn_3d.pack(side=tk.LEFT, padx=3)

# ********************************** FIN DE BOTONES **********************************

# Frecuencia de Excitación
tk.Label(right_panel, text="Frecuencia de Excitación", bg='white', font=('Arial', 9)).pack(pady=(15,2))
label_frecuencia = tk.Label(right_panel, textvariable=Frecuencia, bg='white', 
                           relief=tk.RIDGE, borderwidth=2, width=20, 
                           font=('Arial', 10), anchor='center')
label_frecuencia.pack(padx=10, pady=2)

# Dato 1 para mostrar
tk.Label(right_panel, text="Dato 1 para mostrar", bg='white', font=('Arial', 9)).pack(pady=(10,2))
label_dato1 = tk.Label(right_panel, textvariable=Dato1, bg='white', 
                      relief=tk.RIDGE, borderwidth=2, width=20, 
                      font=('Arial', 10), anchor='center')
label_dato1.pack(padx=10, pady=2)

# Dato 2 para mostrar
tk.Label(right_panel, text="Dato 2 para mostrar", bg='white', font=('Arial', 9)).pack(pady=(10,2))
label_dato2 = tk.Label(right_panel, textvariable=Dato2, bg='white', 
                      relief=tk.RIDGE, borderwidth=2, width=20, 
                      font=('Arial', 10), anchor='center')
label_dato2.pack(padx=10, pady=2)

# Status
tk.Label(right_panel, text="Status", bg='white', font=('Arial', 9)).pack(pady=(10,2))
label_dato3 = tk.Label(right_panel, textvariable=Dato3, bg='white', 
                      relief=tk.RIDGE, borderwidth=2, width=20, 
                      font=('Arial', 10), anchor='center')
label_dato3.pack(padx=10, pady=2)

# ********************************** FIN DE ESTRUCTURA PRINCIPAL **********************************

# Al cerrar la ventana
def on_closing():
    plt.close('all')
    root.destroy()

root.protocol("WM_DELETE_WINDOW", on_closing)

root.mainloop()