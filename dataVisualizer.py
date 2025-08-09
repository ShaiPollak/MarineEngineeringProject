import webbrowser
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
from ncFileReader import RightSideData, LeftSideData

class DataVisualizer(tk.Tk):
    def __init__(self, nc_file_path, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.title('SWOT Project Data Visualizer')
        self.geometry('1200x600')

        self.Z = None
        self.X = None
        self.Y = None

        self.nc_file_path = nc_file_path
        self.side = tk.StringVar()
        self.row = tk.IntVar()
        self.pixel = tk.IntVar()
        self.angle = tk.IntVar()
        self.scale = tk.StringVar(value='Pixels')  # Default scale
        self.colormap_type = tk.StringVar(value='continuous')

        '''
        # Dropdown to select 'left' or 'right'
        side_label = ttk.Label(self, text="Select Side:")
        side_label.pack(pady=10)
        side_dropdown = ttk.Combobox(self, textvariable=self.side, values=['left', 'right'])
        side_dropdown.pack()
        '''

        # Entry to input the row number
        row_label = ttk.Label(self, text="Enter Row Number:")
        row_label.pack(pady=10)
        row_entry = ttk.Entry(self, textvariable=self.row)
        row_entry.pack()

        # Button to open Google Maps
        maps_button = ttk.Button(self, text="Open Google Maps", command=self.open_google_maps)
        maps_button.pack(pady=5)

        # Dropdown for colormap type selection
        colormap_label = ttk.Label(self, text="Select Colormap Type:")
        colormap_label.pack(pady=10)
        colormap_dropdown = ttk.Combobox(self, textvariable=self.colormap_type, values=['continuous', 'discrete'])
        colormap_dropdown.pack()

        # Dropdown for scale selection
        scale_label = ttk.Label(self, text="Select Scale:")
        scale_label.pack(pady=10)
        scale_dropdown = ttk.Combobox(self, textvariable=self.scale, values=['Pixels', 'Meters'])
        scale_dropdown.pack()

        # Button to plot the graph
        plot_button = ttk.Button(self, text="Plot 3D SSH", command=self.plot_ssh)
        plot_button.pack(pady=20)

        # Button to align view to the YZ plane
        yz_button = ttk.Button(self, text="View YZ Plane", command=lambda: self.set_view(0, 90))
        yz_button.pack(pady=5)

        # Button to align view to the XZ plane
        xz_button = ttk.Button(self, text="View XZ Plane", command=lambda: self.set_view(0, 0))
        xz_button.pack(pady=5)

        # Button to align view to the XY plane
        xy_button = ttk.Button(self, text="View XY Plane", command=lambda: self.set_view(90, -90))
        xy_button.pack(pady=5)

        # Label for displaying coordinates
        self.coord_label = ttk.Label(self, text="Coordinates: ")
        self.coord_label.pack(pady=10)

        # Placeholder for matplotlib figure
        self.fig = plt.figure(figsize=(14, 6))
        self.ax1 = self.fig.add_subplot(121, projection='3d')
        self.ax2 = self.fig.add_subplot(122, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)

        self.canvas_widget.bind("<MouseWheel>", self.zoom_fun)


    def withdrawSSHData(self, side):

        SSH_data = None
        selected_row = self.row.get()

        if side == 'left':
            left_data = LeftSideData(self.nc_file_path)
            SSH_data = left_data.get_SSH_matrix(selected_row, selected_row+239, 1, 240, True)
            print(SSH_data)
        if side == 'right':
            right_data = RightSideData(self.nc_file_path)
            SSH_data = right_data.get_SSH_matrix(selected_row, selected_row+239, 1, 240, True)

        return SSH_data


    def plot_ssh(self):


        # Create a 3D plot
        Z_left = self.withdrawSSHData('left')
        Z_right = self.withdrawSSHData('right')

        if Z_left is None or Z_right is None:
            messagebox.showerror("Error", "Failed to retrieve SSH data.")
            return

        scale_factor = 250 if self.scale.get() == 'Meters' else 1
        x = np.arange(1, Z_left.shape[1] + 1) * scale_factor
        y = np.arange(1, Z_left.shape[0] + 1) * scale_factor
        X, Y = np.meshgrid(x, y)

        '''
        scale_factor = 250 if self.scale.get() == 'Meters' else 1
        X = np.arange(1, 240, 1) * scale_factor
        Y = np.arange(1, 240, 1) * scale_factor
        '''

        cmap = 'viridis' if self.colormap_type.get() == 'continuous' else 'viridis_r'  # Choose a discrete colormap if needed

        self.ax1.clear()
        self.ax1.plot_surface(X, Y, Z_left, cmap=cmap)
        self.ax1.set_xlabel('X [m]' if scale_factor > 1 else 'X [pixels]')
        self.ax1.set_ylabel('Y [m]' if scale_factor > 1 else 'Y [pixels]')

        self.ax2.clear()
        self.ax2.plot_surface(X, Y, Z_right, cmap=cmap)
        self.ax2.set_xlabel('X [m]' if scale_factor > 1 else 'X [pixels]')
        self.ax2.set_ylabel('Y [m]' if scale_factor > 1 else 'Y [pixels]')

        self.canvas.draw()

    def zoom_fun(self, event):
        # Zoom in or out based on the direction of the mouse wheel
        if event.delta > 0:
            # Zoom in (reduce zlim)
            scale_factor = 0.9
        else:
            # Zoom out (increase zlim)
            scale_factor = 1.1

        x_min, x_max = self.ax1.get_xlim()
        y_min, y_max = self.ax1.get_ylim()
        z_min, z_max = self.ax1.get_zlim()

        self.ax1.set_xlim([x_min * scale_factor, x_max * scale_factor])
        self.ax1.set_ylim([y_min * scale_factor, y_max * scale_factor])
        self.ax1.set_zlim([z_min * scale_factor, z_max * scale_factor])

        self.ax2.set_xlim([x_min * scale_factor, x_max * scale_factor])
        self.ax2.set_ylim([y_min * scale_factor, y_max * scale_factor])
        self.ax2.set_zlim([z_min * scale_factor, z_max * scale_factor])

        self.canvas.draw()

    def set_view(self, elev, azim):
        """Set the view angle for the 3D plot."""
        self.ax1.view_init(elev=elev, azim=azim)
        self.ax2.view_init(elev=elev, azim=azim)
        self.canvas.draw()

    def on_mouse_move(self, event):
        if event.inaxes == self.ax1 and self.Z is not None:
            try:
                # Mapping the event xdata and ydata to your data array indices
                ix = int(np.interp(event.xdata, self.X[0], np.arange(240)))
                iy = int(np.interp(event.ydata, self.Y[:, 0], np.arange(240)))
                print((ix,iy))
                if 1 <= ix < 240 and 1 <= iy < 240:
                    z = self.Z[iy, ix]  # Fetch the Z value using data indices
                    self.coord_label.config(text=f"Coordinates: x={ix}, y={iy}, z={z:.2f}")
                else:
                    self.coord_label.config(text="Coordinates: Out of plot bounds")
            except Exception as e:  # Catch potential errors in index conversion or out-of-bounds
                self.coord_label.config(text="Coordinates: Error accessing data")
        else:
            self.coord_label.config(text="Coordinates: ")

    def open_google_maps(self):
        left_data = LeftSideData(self.nc_file_path)
        latitude = left_data.get_latitude(self.row.get(), 4)  # Replace with your latitude
        longitude = left_data.get_longitude(self.row.get(), 4)  # Replace with your longitude
        # Validate and adjust the latitude and longitude values

        if latitude < -90 or latitude > 90:
            print("Latitude is out of bounds. It should be between -90 and 90.")
            return

        if longitude < -180 or longitude > 180:
            if longitude > 180:
                longitude -= 360
            elif longitude < -180:
                longitude += 360

        # Ensure the longitude is now within the correct range
        if longitude < -180 or longitude > 180:
            print("Longitude is still out of bounds after adjustment.")
            return

        # Google Maps URL with marker
        zoom_level = 11  # Zoom level defined by z parameter in URL
        url = f"https://www.google.com/maps?q={latitude},{longitude}&z={zoom_level}"
        webbrowser.open(url)


