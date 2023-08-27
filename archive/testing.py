#just testing

class Thing():
	def __init__(self):
		self.thing1 = [1]
		self.thing2 = self.thing1
	def change_thing2(self):
		self.thing2 = [self.thing2[0] + 2]
	def self_test(self, onj):
		return self == onj


this1 = Thing()
print(this1.self_test(this1))


