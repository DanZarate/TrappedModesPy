import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from Source.calculos import calcular_canal

frames = None
X = None
Z = None

root = tk.Tk()
root.title("Interfaz")
root.geometry("1400x700")
#root.iconbitmap("Source/Logo3.ico")
root.config(bg='#bababa')

# Variables globales
Dimension = tk.StringVar(value="2D") #Por defecto 2D

# Variables
Longitud = tk.StringVar(value="")
Ancho = tk.StringVar(value="")
tipo_entrada = tk.StringVar(value="Escoger Entrada")

# Frame principal
main_frame = tk.Frame(root, bg='lightgray', relief=tk.RIDGE, borderwidth=2)
main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Título
title = tk.Label(main_frame, text="Simulador de Ondas para Detección de Modos Atrapados", font=('Arial', 12, 'bold'), bg='lightgray')
title.pack(pady=10)

# Frame contenedor
content_frame = tk.Frame(main_frame, bg='white', relief=tk.RIDGE, borderwidth=2)
content_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Panel izquierdo (Controles - Fondo)
left_panel = tk.Frame(content_frame, bg='#bababa', width=180)
left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
left_panel.pack_propagate(False)

# Entrada A
tk.Label(left_panel, text="Longitud del Canal",bg='white', font=('Arial', 9)).pack(pady=(5,0))

In_Longitud = tk.Entry(left_panel, 
                   textvariable=Longitud, 
                   width=15, 
                   relief=tk.RIDGE, 
                   borderwidth=2, 
                   justify='center')
In_Longitud.pack(pady=5, padx=10)

# Entrada B
tk.Label(left_panel, text="Ancho del Canal", bg='white', font=('Arial', 9)).pack(pady=(5,0))
In_Ancho = tk.Entry(left_panel, 
                textvariable=Ancho, 
                width=15, 
                relief=tk.RIDGE, 
                borderwidth=2, 
                justify='center')
In_Ancho.pack(pady=5, padx=10)

# Tipo de entrada (Combobox)
tk.Label(left_panel, text="Tipo de entrada", bg='white', font=('Arial', 9)).pack(pady=(15,0))
combo_tipo = ttk.Combobox(left_panel, 
                          textvariable=tipo_entrada,
                          values=["Gota de agua", "Otra", "Otra"], 
                          state='readonly', 
                          width=15)
combo_tipo.pack(pady=5, padx=10)


# ********************************** FUNCIONES **********************************
def plot_2d(frames):
    x = np.linspace(-2*np.pi, 2*np.pi, 500)

    global ax1, ax2, ax3, ax4
    global im, ani
    
    # G1: Wave propagation
    fig1.clear()    
    ax1 = fig1.add_subplot(111)

    # Imagen inicial
    im = ax1.imshow(
        frames[0].T,
        origin='lower',
        cmap='ocean',
        extent=[-float(Longitud.get()), float(Longitud.get()), -float(Ancho.get()), float(Ancho.get())],
        aspect='auto',
        interpolation='bilinear'
    )

    ax1.set_title("Solution with constant depth h", fontsize=12, fontweight='bold', pad=10)
    ax1.set_xlabel('X (m)', fontsize=10)
    ax1.set_ylabel('Z (m)', fontsize=10)
    ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

    cbar1 = plt.colorbar(im, ax=ax1, shrink=0.8, aspect=15)
    cbar1.set_label('Surface Height (m)', rotation=270, labelpad=20, fontsize=9)

    fig1.subplots_adjust(left=0.15, bottom=0.18)

    def update(frame_idx):
        im.set_array(frames[frame_idx].T)
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
        
    # G2: cos(x)
    fig2.clear()
    ax2 = fig2.add_subplot(111)
    ax2.plot(x, float(Longitud.get())*np.cos(x), 'r-', linewidth=2, label='A·cos(x)')
    ax2.set_title('cos(x)', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlabel('x')
    ax2.set_ylabel('Amplitud')
    ax2.legend(loc='upper right')
    fig2.subplots_adjust(left=0.15, bottom=0.18)
    canvas2.draw()
        
    # G3: tan(x)
    fig3.clear()
    ax3 = fig3.add_subplot(111)
    y_tan = np.tan(float(Ancho.get())*x)
    y_tan[np.abs(y_tan) > 10] = np.nan
    ax3.plot(x, y_tan, 'g-', linewidth=2, label='A·tan(x)')
    ax3.set_ylim(-5, 5)
    ax3.set_title('tan(x)', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlabel('x')
    ax3.set_ylabel('Amplitud')
    ax3.legend(loc='upper right')
    fig3.subplots_adjust(left=0.15, bottom=0.18)
    canvas3.draw()
        
    # G4: x^2
    fig4.clear()
    ax4 = fig4.add_subplot(111)
    ax4.plot(x, float(Ancho.get())*x**2, 'm-', linewidth=2, label='A·x²')
    ax4.set_title('x²', fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlabel('x')
    ax4.set_ylabel('Amplitud')
    ax4.legend(loc='upper right')
    fig4.subplots_adjust(left=0.15, bottom=0.18)
    canvas4.draw()
    
def plot_3d(frames, X, Z):
    """Grafica funciones en 3D"""
    # Crear malla para gráficas 3D
    vmin = np.min(frames)
    vmax = np.max(frames)

    global ax1, surf1, ani

    fig1.clear()
    ax1 = fig1.add_subplot(111, projection='3d')  

    def configure_axis():
        ax1.set_xlabel('X (m)', fontsize=9, labelpad=8)
        ax1.set_ylabel('Z (m)', fontsize=9, labelpad=8)
        ax1.set_zlabel('ψ(x,z)', fontsize=9, labelpad=8)
        ax1.set_title("3D wave propagation", fontsize=11, pad=15, fontweight='bold')
        ax1.set_facecolor('#f0f8ff')
        ax1.grid(True, alpha=0.2)
        ax1.xaxis.pane.fill = False
        ax1.yaxis.pane.fill = False
        ax1.zaxis.pane.fill = False
        ax1.set_zlim(vmin, vmax)
        ax1.set_box_aspect((2.5, 1.2, 0.4))

    # Surface inicial
    surf1 = ax1.plot_surface(
        X, Z, frames[0],
        cmap='ocean',
        alpha=0.85,
        vmin=vmin, vmax=vmax,
        antialiased=False,
        edgecolor='none',
        rcount=25, ccount=25,
        shade=False
    )

    configure_axis()

    #cbar = fig1.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10, pad=0.1)
    #cbar.set_label('ψ (m)', rotation=270, labelpad=15, fontsize=9)

    #fig1.subplots_adjust(left=0.15, bottom=0.18)

    def update(i):
        ax1.clear()
        surf = ax1.plot_surface(
            X, Z, frames[i],
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
        return surf,


    ani = animation.FuncAnimation(
        fig1,
        update,
        frames=len(frames),
        interval=50, 
        blit=False,
        repeat=False
    )
    

    canvas1.draw()   
        
    # G2: cos(x) en 3D
    fig2.clear()
    ax2_3d = fig2.add_subplot(111, projection='3d')
    ax2_3d.set_title('cos(√(x²+y²))', fontsize=10)
    ax2_3d.set_xlabel('x')
    ax2_3d.set_ylabel('y')
    ax2_3d.set_zlabel('Amplitud')

    canvas2.draw()
        
    # G3: tan(x) en 3D (limitado)
    fig3.clear()
    ax3_3d = fig3.add_subplot(111, projection='3d')
    ax3_3d.set_title('sin(x)·cos(y)', fontsize=10)
    ax3_3d.set_xlabel('x')
    ax3_3d.set_ylabel('y')
    ax3_3d.set_zlabel('Amplitud')

    canvas3.draw()
        
    # G4: x^2 en 3D
    fig4.clear()
    ax4_3d = fig4.add_subplot(111, projection='3d')
    ax4_3d.set_title('x²+y²', fontsize=10)
    ax4_3d.set_xlabel('x')
    ax4_3d.set_ylabel('y')
    ax4_3d.set_zlabel('Amplitud')

    canvas4.draw()

def Calcular():      
    global frames
    global X, Z
    # print("Aquí iría el código")
    # Validamos entradas
    try:
        lon = float(Longitud.get())
        ancho = float(Ancho.get())
        print(f"Longitud del canal: {lon}")
        print(f"Ancho del canal: {ancho}")
        print(f"Tipo de entrada: {tipo_entrada.get()}")

        frames, X, Z = calcular_canal(lon, ancho)
        print("Cálculo completado.")
    except ValueError:
        print("Error: Por favor ingrese valores numéricos válidos para Longitud y Ancho.")
    
def Graficar():
    # Grafica las funciones en 2D o 3D
    if str(Dimension.get()) == "2D":
        plot_2d(frames)
    else:
        plot_3d(frames, X, Z)

def set_dimension(mode):
    # Cambia entre modo 2D y 3D
    Dimension.set(mode)
        
    if str(Dimension.get()) == "2D":
        btn_2d.config(relief=tk.SUNKEN, bg='#20AD4E', borderwidth=3)
        btn_3d.config(relief=tk.RAISED, bg='lightgray', borderwidth=2)
        print("Modo 2D activado")
    else:
        btn_3d.config(relief=tk.SUNKEN, bg='#20AD4E', borderwidth=3)
        btn_2d.config(relief=tk.RAISED, bg='lightgray', borderwidth=2)
        print("Modo 3D activado")

# ********************************** FIN DE FUNCIONES **********************************


# Botón Calcular
btn_calcular = tk.Button(left_panel, text="Calcular", bg="#58DB83", 
    relief=tk.RAISED, borderwidth=2, width=13,
    command= Calcular, font=('Arial', 9, 'bold'))
btn_calcular.pack(pady=10, padx=10)
        
# Botón Graficar
btn_graficar = tk.Button(left_panel, text="Graficar", bg="#58DB83", 
    relief=tk.RAISED, borderwidth=2, width=13,
    command= Graficar, font=('Arial', 9, 'bold'))
btn_graficar.pack(pady=5, padx=10)

# Botones 2D y 3D
dimension_frame = tk.Frame(left_panel, bg='white')
dimension_frame.pack(pady=20, padx=10)
        
# Botón 2D        
btn_2d = tk.Button(dimension_frame, text="2D", bg='#20AD4E', 
    relief=tk.SUNKEN, borderwidth=3, width=6,
    command=lambda: set_dimension("2D"),
    font=('Arial', 10, 'bold'))
btn_2d.pack(pady=5)
        
# Botón 3D
btn_3d = tk.Button(dimension_frame, text="3D", bg='lightgray', 
    relief=tk.RAISED, borderwidth=2, width=6,
    command=lambda: set_dimension("3D"),
    font=('Arial', 10, 'bold'))
btn_3d.pack(pady=5)


# ********************************** ESTRUCTURA GRÁFICAS **********************************

# SubPanel (gráficas)
Subpanel = tk.Frame(content_frame, bg='white')
Subpanel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=10)

# Proporciones SubPanel
Subpanel.rowconfigure(0, weight=1)
Subpanel.columnconfigure(0, weight=1)
Subpanel.columnconfigure(1, weight=1)

# Panel Central (gráficas G1 y G2)
Central_panel = tk.Frame(Subpanel, bg='white')
Central_panel.grid(row=0, column=0, sticky="nsew", padx=5)  

# Proporciones G1 y G2
Central_panel.rowconfigure(0, weight=1)
Central_panel.rowconfigure(1, weight=1)
Central_panel.columnconfigure(0, weight=1)

# Frame para G1
g1_frame = tk.Frame(Central_panel, relief=tk.RIDGE, borderwidth=2)
g1_frame.grid(row=0, column=0, sticky="nsew", pady=(0,5))
tk.Label(g1_frame, text="G1", bg='white', font=('Arial', 10)).pack()
fig1 = Figure(figsize=(4, 2.5), facecolor='white')
ax1 = fig1.add_subplot(111)
canvas1 = FigureCanvasTkAgg(fig1, g1_frame)
canvas1.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# Frame para G2
g2_frame = tk.Frame(Central_panel, relief=tk.RIDGE, borderwidth=2)
g2_frame.grid(row=1, column=0, sticky="nsew")
tk.Label(g2_frame, text="G2", bg='white', font=('Arial', 10)).pack()
fig2 = Figure(figsize=(4, 2.5), facecolor='white')
ax2 = fig2.add_subplot(111)
canvas2 = FigureCanvasTkAgg(fig2, g2_frame)
canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
# Panel derecho (gráficas G3 y G4)
right_panel = tk.Frame(Subpanel, bg='white')
right_panel.grid(row=0, column=1, sticky="nsew", padx=5)

# Proporciones G3 y G4
right_panel.rowconfigure(0, weight=1)
right_panel.rowconfigure(1, weight=1)
right_panel.columnconfigure(0, weight=1)
        
# Gráfica G3
g3_frame = tk.Frame(right_panel, relief=tk.RIDGE, borderwidth=2)
g3_frame.grid(row=0, column=0, sticky="nsew", pady=(0,5))
tk.Label(g3_frame, text="G3", bg='white', font=('Arial', 10)).pack()

fig3 = Figure(figsize=(4, 2.5), facecolor='white')
ax3 = fig3.add_subplot(111)
canvas3 = FigureCanvasTkAgg(fig3, g3_frame)
canvas3.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
# Gráfica G4
g4_frame = tk.Frame(right_panel, relief=tk.RIDGE, borderwidth=2)
g4_frame.grid(row=1, column=0, sticky="nsew")
tk.Label(g4_frame, text="G4", bg='white', font=('Arial', 10)).pack()
fig4 = Figure(figsize=(4, 2.5), facecolor='white')
ax4 = fig4.add_subplot(111)
canvas4 = FigureCanvasTkAgg(fig4, g4_frame)
canvas4.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# ********************************** FIN DE ESTRUCTURA GRÁFICAS **********************************

root.mainloop()