import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg
import matplotlib
matplotlib.use('TkAgg')

import Ts_cycle

"""
Demonstrates one way of embedding Matplotlib figures into a PySimpleGUI window.
Basic steps are:
 * Create a Canvas Element
 * Layout form
 * Display form (NON BLOCKING)
 * Draw plots onto convas
 * Display form (BLOCKING)
 
 Based on information from: https://matplotlib.org/3.1.0/gallery/user_interfaces/embedding_in_tk_sgskip.html
 (Thank you Em-Bo & dirck)
"""

def draw_compressor(graph, pressure_ratio, start_x, axis_y, scale, length):
    # pressure_ratio relates to scaling of compressor
    # start_x is the point on x axis where drawing starts
    # axis_y is the centreline of the compressor
    # scale is the height of the input in px
    length = int(scale*length)
    exit_height = int(0.8 * scale * np.power(pressure_ratio, -0.25))
    points = []
    points.append([start_x, axis_y-(scale/2)])
    points.append([start_x+length, axis_y-(exit_height/2)])
    points.append([start_x+length, axis_y+(exit_height/2)])
    points.append([start_x, axis_y+(scale/2)])
    
    graph.draw_polygon(points, 
        fill_color="#999999", line_color="#000000", line_width=1)

def draw_HX(graph, length, start_x, axis_y, scale):
    draw_length = scale*length

    # External box
    graph.draw_rectangle([start_x, axis_y+(scale/2)],[start_x+draw_length, axis_y-(scale/2)],
        fill_color="#999999", line_color="#000000", line_width=1)

    # Through-line
    graph.draw_line([start_x, axis_y], [start_x+draw_length, axis_y], width=2)

    # Draw zigzags
    zigzag_length_total = draw_length - (scale * 0.2)
    zigzag_start = start_x + (scale*0.1)
    num_zags = int(np.floor(zigzag_length_total / (scale/4)))
    zigzag_length = zigzag_length_total / num_zags

    zigzag_points = []
    zigzag_points.append([zigzag_start, axis_y-(0.75*scale)])
    for n in range(num_zags):
        zigzag_points.append([zigzag_start + (n)*zigzag_length, axis_y+(0.25*scale)])
        zigzag_points.append([zigzag_start + (n+0.5)*zigzag_length, axis_y-(0.25*scale)])

    zigzag_points.append([zigzag_start+zigzag_length_total, axis_y+(0.25*scale)])
    zigzag_points.append([zigzag_start+zigzag_length_total, axis_y-(0.75*scale)])

    for i in range(len(zigzag_points)-1):
        graph.draw_line(zigzag_points[i], zigzag_points[i+1], color="blue")

def draw_compressor_string(graph, machine_params):
    components = []
    graph.erase()
    for i in machine_params:
        comp = []
        comp.append(i["type"])
        if i["type"] == "c":
            comp.append(i["pressure ratio"])
            comp.append(50)
            comp.append(2**(-2.5 + i["pressure ratio"]))
            components.append(comp)
        elif i["type"] == "h":
            comp.append(i["delta h"]/30e3)
            comp.append(50)
            components.append(comp)
            
    #components = [
    #    ["c", 12, 75, 0.5] ,["h", 2, 50], ["c", 12, 75, 0.5] ,["h", 2, 50], ["c", 12, 75, 0.5]
    #]

    start_axis_1 = 50
    axis_offset = 100
    global_start_x = 50

    start_x = 50
    axis_y = start_axis_1

    for c in components:
        if c[0] == "c":
            draw_compressor(graph, c[1], start_x, axis_y, c[2], c[3])
            start_x += c[2]*c[3]
            line_length = (c[2]*c[3]*0.4)
        elif c[0] == "h":
            draw_HX(graph, c[1], start_x, axis_y, c[2])
            start_x += c[1]*c[2]
            line_length = c[1] * c[2] * 0.25
        
        if start_x < 700:
            graph.draw_line([start_x, axis_y], [start_x+line_length, axis_y], width=2)
            start_x += line_length
        else:
            graph.draw_line([start_x, axis_y], [start_x+20, axis_y], width=2)
            graph.draw_line([start_x+20, axis_y], [start_x+20, axis_y+(0.5*axis_offset)], width=2)
            graph.draw_line([start_x+20, axis_y+(0.5*axis_offset)], [global_start_x-20, axis_y+(0.5*axis_offset)], width=2)
            graph.draw_line([global_start_x-20, axis_y+(0.5*axis_offset)], [global_start_x-20, axis_y+axis_offset], width=2)
            graph.draw_line([global_start_x-20, axis_y+axis_offset], [global_start_x, axis_y+axis_offset], width=2)
            
            start_x = global_start_x
            axis_y += axis_offset

# ------------------------------- END OF YOUR MATPLOTLIB CODE -------------------------------

# ------------------------------- Beginning of Matplotlib helper code -----------------------

def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas.TKCanvas)
    canvas.set_size((900,300))

    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg

def update_figure(canvas, figure_canvas_agg):
    figure_canvas_agg.get_tk_widget().forget()
    figure_canvas_agg = draw_figure(canvas, fig)
    return(figure_canvas_agg)

# ------------------------------- Beginning of GUI CODE -------------------------------


visualisation_column = [
    [sg.Canvas(key='-CANVAS-', size=(600, 500), expand_x=True, expand_y=True)],
    [sg.HSeparator()],
    [sg.Graph(canvas_size=(900,300), graph_bottom_left=(0, 0), graph_top_right=(900,300), background_color="#ffffff", key="-DRAW-")]
]

system_performance_column = [
    [sg.Text("Cycle inputs")],
    [sg.Radio("Radiator-driven cycle", "-MODE-SELECT-", key="-MODE-RAD-", enable_events=True, default=True),
    sg.Radio("Icemelt-driven cycle", "-MODE-SELECT-", key="-MODE-ICE-", enable_events=True)],
    [sg.Text("Number of intercoolers"), sg.Input("10", key="-N-ICs-", size=(3,1), enable_events=False)],
    [sg.Text("Input pressure"), sg.Input("900", key="-P-IN-", size=(4,1), enable_events=False), sg.Text("Pa")],
    [sg.Text("Input temperature"), sg.Input("190", key="-T-IN-", size=(4,1), enable_events=False), sg.Text("K")],
    [sg.Text("Output pressure"), sg.Input("6", key="-P-OUT-", size=(4,1), enable_events=False), sg.Text("bar")],
    [sg.Text("Intercooler outlet \n temperature", key="-IC-TEMP-NAME-"), 
    sg.Input("250", key="-IC-TEMP-", size=(4,1), enable_events=False), sg.Text("K")],
    [sg.Text("Intercooler pressure \n drop per K"), sg.Input("0.05", size=(5,1), enable_events=False, key="-IC-P-DROP-"), sg.Text("%/K")],
    [sg.Button("Plot cycle", key="-CYCLE-DRAW-",enable_events=True)],
    [sg.HSeparator()],
    [sg.Text("System performance summary")],
    [sg.Listbox(values=[], enable_events=False, disabled=True, size=(50,30), key="-PERFORMANCE-DATA-")],
    
]

"""layout = [[sg.Text('Cycle Visualisation System')],
        [sg.Canvas(key='-CANVAS-', size=(600, 500), expand_x=True, expand_y=True), sg.VSeparator(), sg.Column(system_performance_column)],
        [sg.HSeparator()],
        [sg.Graph(canvas_size=(1200,250), graph_bottom_left=(0, 0), graph_top_right=(900,300), background_color="#ffffff", key="-DRAW-")]]"""
layout = [
    [sg.Text('Cycle Visualisation System')],
    [sg.Column(visualisation_column), sg.VSeparator(), sg.Column(system_performance_column)]
]


fig,ax,Trange, srange, hgrid, pgrid, dgrid = Ts_cycle.main_load()

# create the form and show it without the plot
window = sg.Window('TurboMars Cycle Visualisation', layout, finalize=True, element_justification='left', resizable=False)
window.Finalize()
#window.Maximize()
sg.set_options(scaling=100 / 72)

#draw_compressor_string(window["-DRAW-"])

fig_canvas_agg = draw_figure(window['-CANVAS-'], fig)

while True:  # Event Loop
    event, values = window.read()
    print(event, values)
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    if event == 'Show':
        # Update the "output" text element to be the value of "input" element
        window['-OUTPUT-'].update(values['-IN-'])
    
    elif event == "-MODE-ICE-":
        window["-IC-TEMP-NAME-"].update(value="Temperature margin \n above sublimation")
        window["-IC-TEMP-"].update(value=40)
    elif event == "-MODE-RAD-":
        window["-IC-TEMP-NAME-"].update(value="Intercooler outlet \n temperature")
        window["-IC-TEMP-"].update(value=250)
        
    
    elif event == "-CYCLE-DRAW-":
        plt.clf()
        plt.cla()
        
        fig, ax = None, None
        fig, ax = Ts_cycle.main_draw(Trange, srange, hgrid, pgrid, dgrid)
        fig.set_size_inches(6, 4)

        try:
            # Load in numeric values from the input fields
            n_ics = int(values["-N-ICs-"])
            p_in = float(values["-P-IN-"])
            T_in = float(values["-T-IN-"])
            p_out = float(values["-P-OUT-"]) * 1e5
            T_IC = float(values["-IC-TEMP-"])
            IC_pdrop = float(values["-IC-P-DROP-"]) / 100

            # Generate the cycle information 
            if values["-MODE-RAD-"]:
                # Radiator-driven cycle
                pressure_ratios, efficiencies, IC_temps, IC_pressure_drops = Ts_cycle.generate_radiator_cycle(n_ics, T_in, p_in, p_out, T_IC, IC_pdrop)

                print("rad")
            elif values["-MODE-ICE-"]:
                # Icemelt-driven cycle
                pressure_ratios, efficiencies, IC_temps, IC_pressure_drops = Ts_cycle.generate_icemelt_cycle(n_ics, T_in, p_in, p_out, T_IC, IC_pdrop)
            
            # Determine all the properties in the cycle based on the actual gas
            machine_params = Ts_cycle.plot_ideal_process(ax, 
                [T_in, p_in], pressure_ratios, efficiencies, IC_temps, IC_pressure_drops, 
                Trange, srange, hgrid, pgrid, dgrid, return_machine_params=True)

            # Draw the compressor visualisation
            draw_compressor_string(window["-DRAW-"], machine_params)
            
            # Write all the compressor information to the Listbox
            comp_count = 1
            hx_count = 1
            properties = []

            # Get master performance properties
            if machine_params[0]["type"] == "cycle":
                properties.append("Total work {:.2f}kJ/kg, total cooling {:.2f}kJ/kg".format(
                    machine_params[0]["total work"]/1000, machine_params[0]["total cooling"]/1000))
                properties.append("Output pressure {:.2f}bar, outlet temperature {:.0f}K".format(
                    machine_params[0]["p out"]/1e5, machine_params[0]["T out"]))
                properties.append("")

            for machine in machine_params:
                if machine["type"] == "c":
                    properties.append("Compressor Block #"+str(comp_count))
                    properties.append("  Inlet pressure {:.3f}kPa, inlet temp {:.0f}K, PR {:.3f}".format(
                        machine["p in"]/1000, machine["T in"], machine["pressure ratio"]))
                    comp_count+=1
                elif machine["type"] == "h":
                    properties.append("Heat Exchanger #"+str(hx_count))
                    properties.append("  Inlet pressure {:.3f}kPa, inlet temp {:.0f}K, cooling {:.2f}kJ/kg".format(
                        machine["p in"]/1000, machine["T in"], machine["delta h"]/1000))
                    hx_count += 1
            
            #fig.set_size_inches(400 * 2 / float(100), 400 / float(100))

            window["-PERFORMANCE-DATA-"].update(values=properties, disabled=False)
            fig_canvas_agg = update_figure(window['-CANVAS-'], fig_canvas_agg)
            
            
            

        except:
            print("Error")
        window["-CANVAS-"].set_size((900,300))
        window.refresh()
        

window.close()