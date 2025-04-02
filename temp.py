import tensorflow as tf
import numpy as np

# Generate synthetic dataset (3 classes, 2 features)
np.random.seed(42)
X_train = np.random.randn(5000, 2)  # 5000 samples, 2 features
true_weights = np.array([[2, -1, 1], [-3, 2, -2]])  # 2 features, 3 output classes
true_bias = np.array([-1, 1, 0])  # 3 classes

logits = np.dot(X_train, true_weights) + true_bias
y_train = np.argmax(logits, axis=1)  # Assign class labels (0, 1, or 2)

# Convert to TensorFlow tensors
X_train_tf = tf.constant(X_train, dtype=tf.float32)
y_train_tf = tf.constant(y_train, dtype=tf.int32)

# Define model parameters (weights and bias)
num_features = X_train.shape[1]
num_classes = 3  # Multi-class classification

weights = tf.Variable(tf.random.normal(shape=(num_features, num_classes)), dtype=tf.float32)
bias = tf.Variable(tf.zeros(shape=(num_classes,)), dtype=tf.float32)

# Softmax Logistic Regression Model
def softmax_regression(X):
	logits = tf.matmul(X, weights) + bias
	return tf.nn.softmax(logits)

# Loss function: Categorical Cross-Entropy
def loss_fn(y_true, y_pred):
	y_true_one_hot = tf.one_hot(y_true, depth=num_classes)  # Convert labels to one-hot
	return tf.reduce_mean(tf.keras.losses.categorical_crossentropy(y_true_one_hot, y_pred))

# Training Parameters
learning_rate = 0.2
epochs = 5000
optimizer = tf.optimizers.SGD(learning_rate)

# Training Loop
for epoch in range(epochs):
	with tf.GradientTape() as tape:
		y_pred = softmax_regression(X_train_tf)
		loss = loss_fn(y_train_tf, y_pred)
	
	# Compute gradients
	gradients = tape.gradient(loss, [weights, bias])
	
	# Update weights using gradient descent
	optimizer.apply_gradients(zip(gradients, [weights, bias]))
	
	if epoch % 50 == 0:
		print(f"Epoch {epoch}, Loss: {loss.numpy():.4f}")

# Final Weights & Bias
print("Learned Weights:\n", weights.numpy())
print("Learned Bias:\n", bias.numpy())