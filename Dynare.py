import numpy as np
import sympy as sy
import statsmodels.api as sm
import matplotlib.pyplot as plt

from numpy.linalg import eig, inv

sy.init_printing(use_latex='mathjax')
np.set_printoptions(precision=4, suppress=True)

from matplotlib import rcParams

# Restore old behavior of rounding default axis ranges
rcParams['axes.autolimit_mode'] = 'round_numbers'
rcParams['axes.xmargin'] = 0
rcParams['axes.ymargin'] = 0

# Adjust tick placement
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True

# Disable legend frame
rcParams['legend.frameon'] = False

class Dynare(object):
	
	def __init__(self, var, varexo, param_values, model, initval):
		self.var = var
		self.varexo = varexo
		self.param_values = param_values
		self.model = model
		self.initval = initval
		
		self.TranslateInputs()
		
		return None
	
	def TranslateInputs(self):
		"""
		Converts input strings into sympy objects
		and translates timings
		"""
		
		# Endogenous variables
		self.var_symbols = sy.symbols(self.var)
		
		self.n = len(self.var_symbols)

		self.current = ()
		self.future  = ()
		self.past    = ()
		self.steadys = ()

		for i in range(self.n):
			self.current += (sy.symbols(str(self.var_symbols[i])+str('_{t}'))),
			self.future  += (sy.symbols(str(self.var_symbols[i])+str('_{t+1}'))),
			self.past    += (sy.symbols(str(self.var_symbols[i])+str('_{t-1}'))),
			self.steadys += (sy.symbols(str(self.var_symbols[i])+str('_{ss}'))),
			
		
		# Exogenous variables
		self.varexo_symbols = sy.symbols(self.varexo)
		
		try:
			self.p = len(self.varexo_symbols)
		except:
			self.varexo_symbols = sy.symbols(self.varexo),
			self.p = len(self.varexo_symbols)
		
		self.shocks = ()
		
		for i in range(self.p):
			self.shocks += (sy.symbols(str(self.varexo_symbols[i])+str('_{t}'))),
			
		# Model equations
		self.symbol_dict = {'betta':sy.symbols('beta'), 'gama':sy.symbols('gamma')}
		self.timing_dict = {}
		
		for i in range(self.n):
			self.timing_dict[sy.sympify(str(self.var_symbols[i])+'(+1)')] = self.future[i]
			self.timing_dict[sy.sympify(str(self.var_symbols[i])+'(-1)')] = self.past[i]
			self.timing_dict[sy.sympify(str(self.var_symbols[i]))]        = self.current[i]
		
		for i in range(self.p):
			self.timing_dict[sy.sympify(str(self.varexo_symbols[i]))] = self.shocks[i]
		
		self.model_symbols = sy.sympify(self.model)
		self.model_symbols = self.model_symbols.subs(self.timing_dict)
		self.model_symbols = self.model_symbols.subs(self.symbol_dict)
		
		try:
			temp = len(self.model_symbols)
		except:
			self.model_symbols = self.model_symbols,
		
		self.system = sy.Matrix(self.model_symbols)
		
		return None
	
	def SteadySystem(self):
		"""
		Removes lead/lag structure from the model to prepare for steady state calculation
		"""
		
		self.steady_vars = {}
		for i in range(self.n):
			self.steady_vars[self.current[i]] = self.steadys[i]
			self.steady_vars[self.future[i]]  = self.steadys[i]
			self.steady_vars[self.past[i]]    = self.steadys[i]
		for i in range(self.p):
			self.steady_vars[self.shocks[i]]  = 0
		
		self.steady_system = self.system.subs(self.steady_vars)
		
		return self.steady_system
	
	def SteadyValues(self):
		"""
		Numerically solves for the steady state of the system given initval
		"""
		
		self.SteadySystem()
		
		try:
			ss = sy.nsolve(self.steady_system.subs(self.param_values), self.steadys, self.initval)
		except:
			raise RuntimeError('Adjust initial values')
		ss = ss.T.tolist()

		self.steady_values = {}
		for i in range(self.n):
			self.steady_values[self.steadys[i]] = np.float(ss[0][i])
				
		return self.steady_values
	
	def steady(self):
		"""
		Prints out steady state values for the user
		"""
		
		self.SteadyValues()
		
		print('\n' + 'STEADY-STATE RESULTS' + '\n')

		for i in range(self.n):
			print(str(self.var_symbols[i]), '\t%.4f' % self.steady_values[self.steadys[i]])
		
		return None
	
	def resid(self):
		
		self.SteadyValues()
		
		temp = self.steady_system.subs(self.param_values).subs(self.steady_values)
		
		print('\n' + 'Residuals of the static equations' + '\n')
		
		for i in range(self.n):
			print('Equation number', i, ': %.4f' % temp[i])
		
		return None
		
	
	def TimeIteration(self):
		"""
		Solves the first-order approximation of the model using time iteration
		(thanks to Pontus Rendahl for teaching me the method)
		"""
		
		self.SteadyValues()
		
		self.A_symb = self.system.jacobian(self.past)
		self.B_symb = self.system.jacobian(self.current)
		self.C_symb = self.system.jacobian(self.future)
		self.D_symb = self.system.jacobian(self.shocks)
		
		self.A = np.array(self.A_symb.subs(self.steady_vars).subs(self.steady_values).subs(self.param_values)).astype(float)
		self.B = np.array(self.B_symb.subs(self.steady_vars).subs(self.steady_values).subs(self.param_values)).astype(float)
		self.C = np.array(self.C_symb.subs(self.steady_vars).subs(self.steady_values).subs(self.param_values)).astype(float)
		self.D = np.array(self.D_symb.subs(self.steady_vars).subs(self.steady_values).subs(self.param_values)).astype(float)
		
		self.metric = 1
		self.F = np.zeros((self.n, self.n))
		self.S = np.zeros((self.n, self.n))

		# Add maxit to while loop?
		
		while self.metric > 1e-13:
			self.F = inv(self.B + self.C @ self.F) @ (-self.A)
			self.S = inv(self.B + self.A @ self.S) @ (-self.C)

			self.metric = np.max(np.max(np.abs(self.A + self.B @ self.F + self.C @ self.F @ self.F)))

		self.Q = -inv(self.B + self.C @ self.F) @ self.D
		
		# Need formal BK check?

		if sum(eig(self.F)[0] > 1) != 0:
			raise RuntimeError('Blanchard Kahn conditions are not satisfied: no stable equilibrium')
		if sum(eig(self.S)[0] > 1) != 0:
			raise RuntimeError('Blanchard Kahn conditions are not satisfied: indeterminacy')
		
		return None
	
	def SimulatedMoments(self, hp_filter=None, shocks_stderr=0.01, periods=10000):
		"""
		Calculates simulated moments
		"""
		
		self.TimeIteration()
		
		x = np.zeros((self.n, periods))
		ɛ = np.zeros((self.p, periods))
		
		for i in range(self.p):
			ɛ[i, :] = shocks_stderr * np.random.randn(periods)

		for t in range(1, periods):
			x[:, t] = self.F @ x[:, t-1] + self.Q @ ɛ[:, t]
		
		print('SIMULATED MOMENTS')
		print('')
			
		if hp_filter == None:
			print('VARIABLE \t STD. DEV.')
			
			for i in range(self.n):
				print(str(self.var_symbols[i]), '\t\t {:.4f}'.format(np.std(x[i, :])))
		else:
			print('VARIABLE \t STD. DEV.')
			
			self.SteadyValues()
			
			hp = np.zeros((self.n, periods))
			try:
				for i in range(self.n):
					hp[i, :], hp_trend = sm.tsa.filters.hpfilter(100*np.log(x[i, :] + self.steady_values[self.steadys[i]]), lamb=hp_filter)
			except:
				print('Error: hp_filter takes only numbers as parameters')
				
			for i in range(self.n):
				print(str(self.var_symbols[i]), '\t\t {:.4f}'.format(np.std(hp[i, :])))
				
			print('')
			print('COEFFICIENTS OF AUTOCORRELATION')
			for i in range(self.n):
				print(str(self.var_symbols[i]), '\t\t {:.4f}'.format(np.corrcoef(hp[i, :-1], hp[i, 1:])[1][0]))
				
			print('')
			print('MATRIX OF CORRELATIONS')
			print('Variables', '\t', str(self.var_symbols[0]))
			for i in range(self.n):
				print(str(self.var_symbols[i]), '\t\t {:.4f}'.format(np.corrcoef(hp[0, :], hp[i, :])[1][0]))
			
			
		
			
		return None
	
	def stoch_simul(self, irf=40, shocks_stderr=None, periods=None):
		"""
		Prints policy and transition functions
		and plots Impluse Response Functions
		"""
		
		self.TimeIteration()

		FT = self.F.T
		QT = self.Q.T

		print('\n'+'POLICY AND TRANSITION FUNCTIONS'+'\n')

		header = '\t'
		for v in self.var_symbols:
			header += '\t' + str(v)
		print(header)
		
		line = ''
		for i in range(self.n):
			line += '\t%.4f' % self.steady_values[self.steadys[i]]
		print('Constant' + line)

		for i in range(self.n):
			if (FT[i] != np.zeros((1, self.n))).any():
				line = '\t'
				for j in range(self.n):
					line += '\t%.4f' % FT[i, j]
				print(str(self.var_symbols[i])+'(-1)', line)
		for i in range(self.p):
			line = '\t'
			for j in range(self.n):
				line += '\t%.4f' % QT[i, j]
			print(str(self.varexo_symbols[i]), '   ', line)
		
		# Impulse response functions
		if irf > 0:
			print('\n')
			x = np.zeros((self.n, irf+2))
			
			for j in range(self.p):
				print('\n\tImpulse response functions to '+str(self.varexo_symbols[j]))
				ɛ = np.zeros((self.p, irf+2))
				if shocks_stderr == None:
					ɛ[j, 1] = 1
				else:
					ɛ[j, 1] = shocks_stderr[j]

				for t in range(1, irf+2):
					x[:, t] = self.F @ x[:, t-1] + self.Q @ ɛ[:, t]
					
				y_dim = int(np.ceil(self.n/3))
				fig, axs = plt.subplots(y_dim, 3, figsize=(16, 4*y_dim), sharex=False, sharey=False)

				for i in range(self.n):
					if sum(abs(x[i, 1:].T)) > 0:
						if self.n <= 3:
							ax = axs[i]
						else:
							ax = axs[i//3, i%3]
						ax.plot(x[i, 1:].T, 'k', lw=2)
						ax.hlines(0, 0, irf, 'r')
						# if self.p == 1:
						# ax.title(str(self.var_symbols[i]))
						ax.set_title(str(self.var_symbols[i]))
						# else:
							# plt.title(str(self.varexo_symbols[j]) + ' -> ' + str(self.var_symbols[i]))
					
					# plt.title('Impulse response functions to '+str(self.varexo_symbols[j]))
				plt.show()
			
		if periods == None:
			return None
		else:
		
			self.TimeIteration()
			
			x = np.zeros((self.n, periods))
			ɛ = np.zeros((self.p, periods))
			
			for i in range(self.p):
				ɛ[i, :] = shocks_stderr * np.random.randn(periods)

			for t in range(1, periods):
				x[:, t] = self.F @ x[:, t-1] + self.Q @ ɛ[:, t]
			
			return x