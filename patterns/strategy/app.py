# en.wikipedia.org/wiki/Strategy_pattern

class Strategy:
	def run(self, a, b):
		raise NotImplementedError


class ConcreteStrategyAdd(Strategy):
	def run(self, a, b):
		return a + b

class ConcreteStrategySubtract(Strategy):
	def run(self, a, b):
		return a - b	

class ConcreteStrategySubtract(Strategy):
	def run(self, a, b):
		return a*b


class Context:
	def __init__(self, strategy):
		self.strategy = strategy
	
	def runStrategy(self, a, b):
		return self.strategy.run(a,b)
		
def main():
	context = Context(ConcreteStrategySubtract())
	print context.runStrategy(1,2)
	
if __name__ == '__main__':
	main()
 
		

