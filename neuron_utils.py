from neuron import h
import random
import numpy as np
import re


class NeuronUtils:
    @staticmethod
    def count_section(section_name):
        return sum([ bool(re.match(rf"^{re.escape(section_name)}", section.name())) for section in h.allsec() ])

    @staticmethod
    def get_compartment(compartment, idx):
        if compartment == "Soma":
            return h.soma[idx]
        if compartment == "Apical":
            return h.apical_dendrite[idx]
        if compartment == "Basal":
            return h.dendrite[idx]

        return None

    @staticmethod
    def get_n_compartments(compartment):
        if compartment == "Soma":
            return NeuronUtils.count_section("soma")
        if compartment == "Apical":
            return NeuronUtils.count_section("apical_dendrite")
        if compartment == "Basal":
            return NeuronUtils.count_section("dendrite")

        return None

    @staticmethod
    def generate_synapse_location(compartment_type, max_distance, min_distance):
        distance_to_soma: int = max_distance + 1
        compartment_idx: int = -1
        loc_pos: float = -1.0
        compartment = None

        n_compartments = NeuronUtils.get_n_compartments(compartment_type)

        while distance_to_soma > max_distance or distance_to_soma < min_distance:
            compartment_idx = random.randint(0, n_compartments - 1)
            loc_pos = random.uniform(0, 1)
            compartment = NeuronUtils.get_compartment(compartment_type, compartment_idx)
            distance_to_soma = h.distance(compartment(loc_pos))

        return compartment, compartment_idx, loc_pos, distance_to_soma

    @staticmethod
    def generate_synapse_location_apical(n_apical_dendrites, max_distance, min_distance, max_diameter):
        """
        Generates a random synapse location within the allowed distance bounds.

        :param n_apical_dendrites: number of apical dendrites
        :param max_distance: maximum distance to soma
        :param min_distance: minimum distance to soma
        :return: dendrite id, position along the dendrite, and distance to the soma
        """

        while True:
            dendrite_idx = random.randint(0, n_apical_dendrites - 1)
            loc_pos = random.uniform(0, 1)
            dendrite = h.apical_dendrite[dendrite_idx]

            distance_to_soma = h.distance(dendrite(loc_pos))
            diameter = dendrite(loc_pos).diam

            if min_distance < distance_to_soma < max_distance and diameter < max_diameter:
                break

        return dendrite_idx, loc_pos, distance_to_soma, diameter

    @staticmethod
    def create_burst_stimulus():
        stimulus = h.BurstStim3()
        stimulus.interval = 25
        stimulus.noise = 0
        stimulus.burstlen = 100
        stimulus.burstint = 250 - stimulus.burstlen

        return stimulus

    @staticmethod
    def create_constant_freq_stimulus():
        frequency = 50
        noise = 0
        delay = 0

        stimulus = h.NetStim()
        stimulus.interval = int(1000 / frequency)
        stimulus.noise = noise
        stimulus.start = delay

        return stimulus

    @staticmethod
    def create_no_stimulus():
        return None

    @staticmethod
    def create_single_stimulus():
        stimulus = h.NetStim()
        stimulus.interval = 1000
        stimulus.noise = 0
        stimulus.start = 50
        stimulus.number = 1

        return stimulus

    @staticmethod
    def set_stimulus(delay = 100, duration = 900, frequency1 = 1000, frequency2 = 1005, amplitude = 100, phase = 0):
        h.setstim(delay, duration, frequency1, frequency2, amplitude, phase)

    @staticmethod
    def set_field(psi_deg = 90, phi_deg = 90):
        field_length = 100
        h.changefield(field_length, psi_deg, phi_deg)

    @staticmethod
    def calculate_synapse_angles(data):
        normal = np.array([h.sField.x3d(1), h.sField.y3d(1), h.sField.z3d(1)])
        normal = normal / np.linalg.norm(normal)

        synapse_angles = []
        for i, synapse in data.iterrows():
            dendrite = int(synapse.Dendrite)
            location = synapse.Location

            i_location = int(h.apical_dendrite[dendrite].n3d() * location)

            max_location = np.floor(h.apical_dendrite[dendrite].n3d()) - 1
            min_location = 0

            i_start = int(np.max([min_location, i_location - 1]))
            i_end = int(np.min([max_location, i_location + 1]))

            dendrite_vector = np.array([
                h.apical_dendrite[dendrite].x3d(i_end) - h.apical_dendrite[dendrite].x3d(i_start),
                h.apical_dendrite[dendrite].y3d(i_end) - h.apical_dendrite[dendrite].y3d(i_start),
                h.apical_dendrite[dendrite].z3d(i_end) - h.apical_dendrite[dendrite].z3d(i_start)
            ])

            dendrite_vector = dendrite_vector / np.linalg.norm(dendrite_vector)

            dot = np.dot(dendrite_vector, normal)
            synapse_angles.append(np.arccos(np.clip(dot, a_min = -1.0, a_max = 1.0)) * (180 / np.pi))

        return synapse_angles
