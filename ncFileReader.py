import math

from netCDF4 import Dataset
import numpy as np
from scipy.stats import linregress

import matplotlib.pyplot as plt

#General File Info Parameters....
class GeneralFileInfo:
    def __init__(self, nc_file_path):

        #Normilze the Slope?
        #self.fix_slope = bool_fix_slope

        #Other Properties
        self.title = None
        self.cycle_number = None
        self.pass_number = None
        self.time_coverage_start = None
        self.time_coverage_end = None
        self.geospatial_lon_min = None
        self.geospatial_lon_max = None
        self.geospatial_lat_min = None
        self.geospatial_lat_max = None
        self.wavelength = None
        self.equator_time = None
        self.equator_longitude = None


        try:
            # Open the NetCDF file
            with Dataset(nc_file_path, "r") as data:

                self.title = data.title
                self.cycle_number = data.cycle_number
                self.pass_number = data.pass_number
                self.time_coverage_start = data.time_coverage_start
                self.time_coverage_end = data.time_coverage_end
                self.geospatial_lon_min = data.geospatial_lon_min
                self.geospatial_lon_max = data.geospatial_lon_max
                self.geospatial_lat_min = data.geospatial_lat_min
                self.geospatial_lat_max = data.geospatial_lat_max
                self.wavelength = data.wavelength
                self.equator_time = data.equator_time
                self.equator_longitude = data.equator_longitude

        except Exception as e:
            print("An error occurred:", e)
    '''
    def get_title(self):
        return self.__title

    def get_cycle_number(self):
        return self.__cycle_number

    def get_pass_number(self):
        return self.__pass_number

    def get_time_coverage_start(self):
        return self.__time_coverage_start

    def get_time_coverage_end(self):
        return self.__time_coverage_end

    def get_geospatial_lon_min(self):
        return self.__geospatial_lon_min

    def get_geospatial_lon_max(self):
        return self.__geospatial_lon_max

    def get_geospatial_lat_min(self):
        return self.__geospatial_lat_min

    def get_geospatial_lat_max(self):
        return self.__geospatial_lat_max

    def get_wavelength(self):
        return self.__wavelength

    def get_equator_time(self):
        return self.__equator_time

    def get_equator_longitude(self):
        return self.__equator_longitude
    '''

class SideData(GeneralFileInfo):
    def __init__(self, nc_file_path, side):
        super().__init__(nc_file_path)
        self.side = side #Right or Left


        self.description = None
        self.num_of_lines = None
        self.num_of_pixels = None

        self.latitude_units = None
        self.latitude_scale_factor = None
        self.__latitude = None

        self.longitude_units = None
        self.longitude_scale_factor = None
        self.__longitude = None

        self.time_units = None
        self.time_scale_factor = None
        self.__time = None

        self.time_tai_units = None
        self.time_tai_scale_factor = None
        self.__time_tai = None

        self.latitude_uncert_units = None
        self.latitude_uncert_scale_factor = None
        self.__latitude_uncert = None

        self.longitude_uncert_units = None
        self.longitude_uncert_scale_factor = None
        self.__longitude_uncert = None

        self.polarization_karin = None

        self.ssh_karin_2_units = None
        self.ssh_karin_2_scale_factor = None
        self.__ssh_karin_2 = None

        self.ssh_karin_2_qual = None

        self.ssh_karin_uncert_units = None
        self.ssh_karin_uncert_scale_factor = None
        self.__ssh_karin_uncert = None

        self.sig0_karin_2_units = None
        self.__sig0_karin_2 = None

        self.__sig0_karin_2_qual = None

        self.sig0_karin_uncert_units = None
        self.__sig0_karin_uncert = None

        self.total_coherence_units = None
        self.total_coherence_scale_factor = None
        self.__total_coherence = None

        self.mean_sea_surface_cnescls_units = None
        self.mean_sea_surface_cnescls_scale_factor = None
        self.__mean_sea_surface_cnescls = None

        self.miti_power_250m_units = None
        self.__miti_power_250m = None

        self.miti_power_var_250m_units = None
        self.__miti_power_var_250m = None

        self.__ancillary_surface_classification_flag = None


        try:
            # Open the NetCDF file
            with Dataset(nc_file_path, "r") as data:

                group = data.groups[self.side]

                self.description = group.description
                self.num_of_lines = group.dimensions['num_lines'].size
                self.num_of_pixels = group.dimensions['num_pixels'].size

                # Attributes for latitude
                self.latitude_units = group.variables['latitude'].units
                self.latitude_scale_factor = group.variables['latitude'].scale_factor
                self.__latitude = group.variables['latitude'][:, :]

                # Attributes for longitude
                self.longitude_units = group.variables['longitude'].units
                self.longitude_scale_factor = group.variables['longitude'].scale_factor
                self.__longitude = group.variables['longitude'][:, :]

                # Attributes for time
                self.time_units = group.variables['time'].units
                self.time_scale_factor = getattr(group.variables['time'], 'scale_factor',
                                                 1)  # Default to 1 if not present
                self.__time = group.variables['time'][:]

                # Attributes for time_tai
                self.time_tai_units = group.variables['time_tai'].units
                self.time_tai_scale_factor = getattr(group.variables['time_tai'], 'scale_factor', 1)
                self.__time_tai = group.variables['time_tai'][:]

                '''
                ------- NOT IMPORTANT FOR NOW ---------
                
                # Attributes for latitude_uncert
                self.latitude_uncert_units = group.variables['latitude_uncert'].units
                self.latitude_uncert_scale_factor = group.variables['latitude_uncert'].scale_factor
                self.latitude_uncert = group.variables['latitude_uncert'][:, :]

                # Attributes for longitude_uncert
                self.longitude_uncert_units = group.variables['longitude_uncert'].units
                self.longitude_uncert_scale_factor = group.variables['longitude_uncert'].scale_factor
                self.longitude_uncert = group.variables['longitude_uncert'][:, :]

                # Attributes for polarization_karin
                self.polarization_karin = group.variables['polarization_karin'][:]
                '''

                # Attributes for ssh_karin_2
                self.ssh_karin_2_units = group.variables['ssh_karin_2'].units
                self.ssh_karin_2_scale_factor = group.variables['ssh_karin_2'].scale_factor
                self.__ssh_karin_2 = group.variables['ssh_karin_2'][:, :]

                '''
                # Attributes for ssh_karin_2_qual
                self.ssh_karin_2_qual = group.variables['ssh_karin_2_qual'][:, :]
                '''

                # Attributes for ssh_karin_uncert
                self.ssh_karin_uncert_units = group.variables['ssh_karin_uncert'].units
                self.ssh_karin_uncert_scale_factor = group.variables['ssh_karin_uncert'].scale_factor
                self.__ssh_karin_uncert = group.variables['ssh_karin_uncert'][:, :]

                '''
                # Attributes for sig0_karin_2
                self.sig0_karin_2_units = group.variables['sig0_karin_2'].units
                self.__sig0_karin_2 = group.variables['sig0_karin_2'][:, :]
                
                
                # Attributes for sig0_karin_2_qual
                self.__sig0_karin_2_qual = group.variables['sig0_karin_2_qual'][:, :]
                
                
                # Attributes for sig0_karin_uncert
                self.sig0_karin_uncert_units = group.variables['sig0_karin_uncert'].units
                self.__sig0_karin_uncert = group.variables['sig0_karin_uncert'][:, :]
                
                # Attributes for total_coherence
                self.total_coherence_units = group.variables['total_coherence'].units
                self.total_coherence_scale_factor = group.variables['total_coherence'].scale_factor
                self.__total_coherence = group.variables['total_coherence'][:, :]
                                '''


                # Attributes for mean_sea_surface_cnescls
                self.mean_sea_surface_cnescls_units = group.variables['mean_sea_surface_cnescls'].units
                self.mean_sea_surface_cnescls_scale_factor = group.variables[
                    'mean_sea_surface_cnescls'].scale_factor
                self.__mean_sea_surface_cnescls = group.variables['mean_sea_surface_cnescls'][:, :]

                '''
                # Attributes for miti_power_250m
                self.miti_power_250m_units = group.variables['miti_power_250m'].units
                self.miti_power_250m = group.variables['miti_power_250m'][:, :]

                # Attributes for miti_power_var_250m
                self.miti_power_var_250m_units = group.variables['miti_power_var_250m'].units
                self.miti_power_var_250m = group.variables['miti_power_var_250m'][:, :]

                # Attributes for ancillary_surface_classification_flag
                self.ancillary_surface_classification_flag = group.variables[
                                                                 'ancillary_surface_classification_flag'][:, :]
                
                '''

        except Exception as e:
            print("An error occurred:", e)

    def get_SSH_matrix(self, start_line, end_line, start_pixel, end_pixel, fix_slope):
        matrix = self.__ssh_karin_2[start_line:end_line + 1, start_pixel:end_pixel + 1]


        #fix nan values and normalize
        for i in range(0, len(matrix)):

            if fix_slope == True:
                matrix[i] = self.rotate_SSH_array(matrix[i])

            else:
                #matrix[i] = self.fill_in_nan(matrix[i])
                pass

        return matrix

    def get_SSH(self, line, pixel): #NORMALIZED!
        tmp = self.fill_in_nan(self.get_SSH_line(line))
        return tmp[pixel]

    def get_SSH_line(self, line):
        #print(self.__ssh_karin_2[line, :])
        fixed_line = self.rotate_SSH_array(self.__ssh_karin_2[line, :])
        #return self.__ssh_karin_2[line, :]
        return fixed_line


    def get_SSH_array_with_angle(self, starting_line, angle_in_deg):
        '''

        :param starting_line: from which line to define 240x240 square?
        :param angle_in_deg: Angle from the line to the hypotenuse
        :return: an array of SSH in this angle, starting from the left bottom corner pixel
        '''

        angle = math.radians(angle_in_deg)
        matrix = self.get_SSH_matrix(starting_line, starting_line+239, 0, 239, True )
        SSH_array = []

        # Case angle < pi/4 (45 deg)
        if (angle < np.pi / 4):
            for i in range(0, 240):
                # calculate the opposite side of the triangle
                btw_row = np.tan(angle) * (i)

                # Cal ceil and floor to that value
                ceil_row = math.ceil(btw_row)
                flr_row = math.floor(btw_row)

                # Get ceil line value and flr line value (in index i)
                SSH1 = matrix[ceil_row, i]
                SSH2 = matrix[flr_row, i]

                # Make an proportional avg
                ceil_proportion = ceil_row - btw_row
                floor_proportion = btw_row - flr_row

                if ceil_proportion == 0:
                    SSH_array.append(SSH1)
                    #print("propotion 0")
                    continue

                proportional_avgSSH = (ceil_proportion * SSH1 + floor_proportion * SSH2)

                # append to list
                SSH_array.append(proportional_avgSSH)

            # Case angle = pi/4 (45 deg): Make just a simple diagonal of the 240x240 square
        if (angle == np.pi / 4):
            for i in range(0, 240):
                SSH = matrix[i, i]

                #Add to the SSH_array list
                SSH_array.append(SSH)

        # Case  pi/2 (90 deg) >= angle > pi/4 (45 deg)
        if ((angle > np.pi / 4) and (angle <= np.pi/2)):
            for i in range(0, 240):
                # calculate the opposite side of the triangle
                btw_pixel = np.tan(np.pi / 2 - angle) * (i)

                # Cal ceil and floor to that value
                ceil_pixel = math.ceil(btw_pixel)
                flr_pixel = math.floor(btw_pixel)

                # Get ceil pixel value and flr pixel value (in index i)
                SSH1 = matrix[i, ceil_pixel]
                SSH2 = matrix[i, flr_pixel]

                # Make an proportional avg
                ceil_proportion = ceil_pixel - btw_pixel
                floor_proportion = btw_pixel - flr_pixel

                if ceil_proportion == 0:
                    SSH_array.append(SSH1)
                    continue

                proportional_avgSSH = (ceil_proportion * SSH1 + floor_proportion * SSH2)

                # append to list
                SSH_array.append(proportional_avgSSH)


        # case 135deg > angle > 90deg, start from the right pixel
        if np.pi / 2 < angle and angle < ((3* np.pi) /4) :
            for i in range(0, 240):
                # Calculate the opposite side of the triangle, adjusting for angles in the 2nd quadrant
                btw_pixel = np.tan(-(np.pi/2 - angle)) * (i)

                # Calculate ceil and floor to that value
                ceil_pixel = math.ceil(btw_pixel)
                flr_pixel = math.floor(btw_pixel)

                # Calculate indices; adjusting the start from pixel 240
                ceil_pixel_index = 239 - ceil_pixel
                flr_pixel_index = 239 - flr_pixel

                # Get ceil pixel value and flr pixel value (in index i)
                SSH1 = matrix[i, ceil_pixel_index]
                SSH2 = matrix[i, flr_pixel_index]

                # Make a proportional average
                ceil_proportion = ceil_pixel - btw_pixel
                floor_proportion = btw_pixel - flr_pixel

                if ceil_proportion == 0:
                    SSH_array.append(SSH1)
                    continue

                proportional_avgSSH = (ceil_proportion * SSH1 + floor_proportion * SSH2)
                # Append to list
                SSH_array.append(proportional_avgSSH)

        # case angle = 135 start from the right pixel
        if (angle == 3*np.pi / 4):
            for i in range(0, 240):
                SSH = matrix[i, 239 - i]

                #Add SSH to list
                SSH_array.append(SSH)

        # case 180 deg > angle > 135deg, start from the right pixel
        if ((3 * np.pi) / 4) < angle and angle < np.pi:
            for i in range(0, 240):
                # calculate the opposite side of the triangle
                btw_row = np.tan((np.pi - angle)) * (i)

                # Cal ceil and floor to that value
                ceil_row = math.ceil(btw_row)
                flr_row = math.floor(btw_row)

                # Get ceil line value and flr line value (in index i)
                SSH1 = matrix[ceil_row, 239 - i]
                SSH2 = matrix[flr_row, 239 - i]

                # Make an proportional avg
                ceil_proportion = ceil_row - btw_row
                floor_proportion = btw_row - flr_row

                if ceil_proportion == 0:
                    SSH_array.append(SSH1)
                    continue

                proportional_avgSSH = (ceil_proportion * SSH1 + floor_proportion * SSH2)

                # append to list
                SSH_array.append(proportional_avgSSH)

        return SSH_array

    def rotate_SSH_array(self, SSH_array):
        #print(SSH_array)
        # Convert None to NaN
        ssh = np.array([np.nan if x is None else x for x in SSH_array], dtype=np.float64)
        #print(ssh)
        # Generate indices for the SSH array
        indices = np.arange(len(ssh))
        # Mask for valid (non-NaN) values
        mask = ~np.isnan(ssh)
        # Fit a linear model to valid SSH values as a function of valid indices
        slope, intercept, r_value, p_value, std_err = linregress(indices[mask], ssh[mask])
        # Calculate the linear trend for all indices
        linear_trend = slope * indices + intercept
        # Detrend the SSH data, keeping NaNs in place
        ssh_detrended = np.where(mask, ssh - linear_trend, np.nan)

        #print(ssh_detrended)

        #Use Adara fix functions to fill the nan values
        ssh_fixed = self.fill_in_nan(ssh_detrended)

        '''
        #DEBUG#
        # Plot the results for visualization
        plt.figure(figsize=(12, 6))

        plt.subplot(2, 1, 1)
        plt.plot(indices, SSH_array, label='Original SSH')
        plt.plot(indices, linear_trend, label='Linear Trend', linestyle='--')
        plt.xlabel('Index')
        plt.ylabel('SSH')
        plt.legend()

        plt.subplot(2, 1, 2)
        plt.plot(indices, ssh_detrended, label='Detrended SSH')
        plt.xlabel('Index')
        plt.ylabel('Detrended SSH')
        plt.legend()

        plt.tight_layout()
        plt.show()
        '''
        #print(ssh_fixed)
        return ssh_fixed

    # Latitude
    def get_latitude_matrix(self, start_line, end_line, start_pixel, end_pixel):
        return self.__latitude[start_line:end_line+1, start_pixel:end_pixel+1]

    def get_latitude(self, line, pixel):
        return self.__latitude[line, pixel]

    def get_latitude_line(self, line):
        return self.__latitude[line, :]

    # Longitude
    def get_longitude_matrix(self, start_line, end_line, start_pixel, end_pixel):
        return self.__longitude[start_line:end_line+1, start_pixel:end_pixel+1]

    def get_longitude(self, line, pixel):
        return self.__longitude[line, pixel]

    def get_longitude_line(self, line):
        return self.__longitude[line, :]


    # Methods for 'time' variable
    def get_time_matrix(self, start_line, end_line):
        return self.__time[start_line:end_line+1]

    def get_time(self, line):
        return self.__time[line]

        # Methods for 'time_tai' variable

    def get_time_tai_matrix(self, start_line, end_line):
        return self.__time_tai[start_line:end_line+1]

    def get_time_tai(self, line):
        return self.__time_tai[line]

    # SSH_Karin Uncertainty
    def get_ssh_karin_uncert_matrix(self, start_line, end_line, start_pixel, end_pixel):
        return self.__ssh_karin_uncert[start_line:end_line+1, start_pixel:end_pixel+1]

    def get_ssh_karin_uncert(self, line, pixel):
        return self.__ssh_karin_uncert[line, pixel]

    def get_ssh_karin_uncert_line(self, line):
        return self.__ssh_karin_uncert[line, :]


    #### OTHER TOOLS #####

    def get_coordinate_of_pixel(self, line, pixel):
        return (self.get_latitude(line, pixel), self.get_longitude(line, pixel))

    def fill_in_nan(self, ssh_array):
        # Fill in NaN values with ssh mean
        mean_ssh = np.nanmean(ssh_array)  # Compute mean ignoring NaNs
        new_ssh = np.where(np.isnan(ssh_array), mean_ssh, ssh_array)  # Fill NaNs with the mean

        # Normalize signal by mean
        new_ssh -= mean_ssh
        return new_ssh


class LeftSideData(SideData):
    def __init__(self, nc_file_path):
        super().__init__(nc_file_path, 'left')

class RightSideData(SideData):
    def __init__(self, nc_file_path):
        super().__init__(nc_file_path, 'right')