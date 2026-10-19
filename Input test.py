import numpy as np
# simulate squiggle 2000 virus data flows
ref = np.random.randn(2000)
test = np.random.randn(2000)
print("ref shape:", ref.shape)
print("test shape:", test.shape)
print("first 5 ref values:", ref[:5])