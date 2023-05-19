import json

f = open(r"C:\Users\couvr\Documents\Justin's Research\TheNextOne\BIOE230_final_work_in_progress.ipynb")

data = json.load(f)

for cell in data['cells']:
	try:
		with open('temp.txt', 'a') as file:
			for line in cell['source']:
				file.write(line)
	except:
		pass