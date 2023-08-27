class temp_Walker():
	def __init__(self):
		self.starting_position = [0,0]
		self.finished = False
		self.positions = []

def load_exploration(raw_data_csv):
	'''import data from an already existing exploration simulation to use in exploitation'''
	#get nutrient level from log file
	nutrient_level = -1
	date_time = raw_data_csv.split('\\')[-2]
	date = date_time.split('_')[0] + '/' + date_time.split('_')[1] + '/' + date_time.split('_')[2]
	time = date_time.split('_')[3].split('-')[0] + ':' + date_time.split('_')[3].split('-')[1]
	with open('log.csv', 'r') as log_file:
		for line in log_file.readlines():
			list_line = line.strip().split(',')
			if (list_line[0] == date) and (time in list_line[1]):
				nutrient_level = int(list_line[3])
				break

	with open(raw_data_csv, 'r') as file:
		header = True
		walker_id = -1
		last_walker_id = -1
		all_walkers = []
		for line in file.readlines():
			if header:
				header = False
				continue

			list_line = line.strip().split(',')
			walker_id = int(list_line[0])
			x = float(list_line[2])
			y = float(list_line[3])

			#starting position
			if walker_id != last_walker_id:
				all_walkers.append(temp_Walker())
				all_walkers[walker_id].starting_position = [x,y]
				all_walkers[walker_id].positions = []
			#positions
			all_walkers[walker_id].positions.append([x,y])

			#finished
			if y > nutrient_level:
				all_walkers[walker_id].finished = True
			else:
				all_walkers[walker_id].finished = False

			last_walker_id = walker_id

	return all_walkers


