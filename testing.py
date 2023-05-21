#just testing

class Thing():
	def __init__(self):
		self.thing1 = [1]
		self.thing2 = self.thing1
	def change_thing2(self):
		self.thing2 = [self.thing2[0] + 2]


list_thing = [Thing() for i in range(2)]

for i in range(2):
	for this in list_thing:
		this.change_thing2()
		print(this.__dict__)


