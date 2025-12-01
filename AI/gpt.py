import numpy as np
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from tcn import TCN
import matplotlib.pyplot as plt


# -------------------------------------------------------------
# 1. Création d'un signal synthétique avec 3 états + bruit
# -------------------------------------------------------------
np.random.seed(0)


segment_lengths = [300, 400, 350]
states = [0, 1, 2]
values = [0.0, 3.0, 1.5]


signal = []
true_states = []


for length, st, val in zip(segment_lengths, states, values):
    segment = np.ones(length) * val + np.random.normal(0, 0.25, length)
    signal.append(segment)
    true_states.append(np.ones(length, dtype=int)*st)


signal = np.concatenate(signal)
true_states = np.concatenate(true_states)


# -------------------------------------------------------------
# 2. Normalisation
# -------------------------------------------------------------
scaler = StandardScaler()
signal_norm = scaler.fit_transform(signal.reshape(-1, 1)).flatten()


# -------------------------------------------------------------
# 3. Construction du dataset (fenêtres temporelles)
# -------------------------------------------------------------
WINDOW = 50


def create_dataset(signal, labels, window):
    X, y = [], []
    for i in range(len(signal) - window):
        X.append(signal[i:i+window])
        y.append(labels[i+window])
        return np.array(X), np.array(y)


X, y = create_dataset(signal_norm, true_states, WINDOW)
X = X[..., np.newaxis] # reshape pour LSTM/TCN


# Séparation train / test
split = int(0.8 * len(X))
X_train, X_test = X[:split], X[split:]
y_train, y_test = y[:split], y[split:]


# -------------------------------------------------------------
# 4. Construction du modèle TCN
# -------------------------------------------------------------
N_STATES = 3


model = tf.keras.Sequential([
tf.keras.layers.Input(shape=(WINDOW, 1)),
TCN(nb_filters=32, kernel_size=3, dilations=[1, 2, 4, 8], dropout_rate=0.1),
tf.keras.layers.Dense(32, activation='relu'),
tf.keras.layers.Dense(N_STATES, activation='softmax')
])


model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
model.summary()
plt.show()
