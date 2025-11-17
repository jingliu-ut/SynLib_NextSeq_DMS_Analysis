import os

# Set the environment variable to use GPU 1
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"

import pandas as pd
import numpy as np
import tensorflow_probability as tfp
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.special import softmax
import argparse
import pickle
import gc
import tensorflow as tf

# Disable GPU
# os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
# os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

tf.keras.backend.set_floatx('float32')

parser = argparse.ArgumentParser('Parameter estimation of FACS-seq')
parser.add_argument('-p', '--pi', help='Propotion of each sample in original library', required=True)
parser.add_argument('-f', '--frac', help='Binned distribution', required=True)
parser.add_argument('-d', '--data', help='cytometry data of the library', required=True)
parser.add_argument('-b', '--boundary', help='FACS boundaries', required=True)
parser.add_argument('-o', '--output', help='output directory', required=True)
parser.add_argument('-c', '--condition', help='condition', required=True)
args = parser.parse_args()
pi = pd.read_csv(args.pi, header=None)
frac = pd.read_csv(args.frac, index_col=0)
data = pd.read_csv(args.data, header=None)
gate = pd.read_csv(args.boundary, header=None)
out_dir = args.output
condition = args.condition

BUFFER_SIZE = 1000000

if condition == "dark" or condition == "light":
    BATCH_SIZE = 3600

elif condition == "ymScarlet":
    BATCH_SIZE = 2400

pi = tf.cast(pi[1].values, tf.float32)
frac = tf.cast(frac.values, tf.float32)
data = tf.cast(data[0].values, tf.float32)
data = tf.reshape(data, [-1, 1])
gate = tf.cast(gate[1].values, tf.float32)
FACS_data = tf.data.Dataset.from_tensor_slices(data).repeat(2).shuffle(BUFFER_SIZE).batch(BATCH_SIZE,
                                                                                          drop_remainder=True)
n_comp = 2

class FS_generator(tf.keras.Model):
    def __init__(self, pi, gate, batch_size):
        super().__init__()
        self.pi = pi
        self.gate = gate
        self.batch_size = batch_size
        mu_mean = None
        sigma_mean = None
        sigma_min = None
        if condition == "dark":
            mu_mean = 2.89
            sigma_mean = 0.26
            sigma_min = 0.09
        elif condition == "light":
            mu_mean = 3.12
            sigma_mean = 0.4
            sigma_min = 0.16
        elif condition == "ymScarlet":
            mu_mean = 3.17
            sigma_mean = 0.92
            sigma_min = 0.32
        # Initialize variables with float32 precision
        self.lamb = tf.Variable(tf.random.normal([pi.shape[0], n_comp], mean=0.5, stddev=0.1, dtype=tf.float32),
                                trainable=True)
        self.mu = tf.Variable(tf.random.normal([pi.shape[0], n_comp], mean=mu_mean, stddev=0.4, dtype=tf.float32),
                              trainable=True)
        self.sigma = tf.Variable(tf.random.normal([pi.shape[0], n_comp], mean=sigma_mean, stddev=0.1, dtype=tf.float32),
                                 trainable=True,
                                 constraint=lambda x: tf.clip_by_value(x, sigma_min, 2.0, 0))

    def __call__(self):
        # Generate samples
        dist = tfp.distributions.MixtureSameFamily(mixture_distribution=tfp.distributions.Categorical(probs=self.pi),
                                                   components_distribution=tfp.distributions.MixtureSameFamily(
                                                       mixture_distribution=tfp.distributions.Categorical(
                                                           logits=self.lamb),
                                                       components_distribution=tfp.distributions.Normal(loc=self.mu,
                                                                                                        scale=self.sigma)))
        sample = dist.sample(sample_shape=[self.batch_size, 1])
        # Calculate cumulative distribution function
        dist_ = tfp.distributions.MixtureSameFamily(
            mixture_distribution=tfp.distributions.Categorical(logits=self.lamb),
            components_distribution=tfp.distributions.Normal(loc=self.mu, scale=self.sigma))
        dist_cdf = tf.concat([[dist_.cdf(g)] for g in self.gate.numpy()], axis=0)
        cdf_diff = tf.concat([tf.reshape(dist_cdf[0], [1, -1]), dist_cdf[1:] - dist_cdf[:-1]], axis=0)
        cdf_diff_tail = tf.reshape(1 - tf.reduce_sum(cdf_diff, axis=0), [1, -1])
        cdf_diff = tf.transpose(tf.concat([cdf_diff, cdf_diff_tail], axis=0))
        return sample, cdf_diff

generator = FS_generator(pi, gate, BATCH_SIZE)

def FS_discriminator():
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(128, activation=tf.nn.leaky_relu))
    model.add(tf.keras.layers.Dense(128, activation=tf.nn.leaky_relu))
    model.add(tf.keras.layers.Dense(1))
    return model

discriminator = FS_discriminator()

# loss function
BinaryCrossentropy = tf.keras.losses.BinaryCrossentropy(from_logits=True)

def discriminator_loss(real_output, fake_output):
    real_loss = BinaryCrossentropy(tf.ones_like(real_output), real_output)
    fake_loss = BinaryCrossentropy(tf.zeros_like(fake_output), fake_output)
    total_loss = real_loss + fake_loss
    return total_loss

def generator_loss(fake_output, cdf_diff):
    dist_loss = tf.math.reduce_sum(pi * tf.keras.losses.categorical_crossentropy(frac, cdf_diff))
    gen_loss = BinaryCrossentropy(tf.ones_like(fake_output), fake_output)
    return dist_loss, gen_loss

generator_optimizer = tf.keras.optimizers.Adam(learning_rate=0.1)
discriminator_optimizer = tf.keras.optimizers.RMSprop(learning_rate=0.01)

def distribution_plot(lamb, mu, sigma, image_name):
    x = np.linspace(-2, 6, 1000)
    fig, ax = plt.subplots()
    ax.hist(np.squeeze(data), bins=x, density=True, color='r', alpha=0.5, label='real_distribution')
    prob = softmax(lamb, axis=1)
    y = [pi.numpy().dot((prob * norm.pdf(i, loc=mu, scale=sigma)).sum(axis=1)) for i in x]
    ax.plot(x, y, color='b', alpha=0.5, label='generated_distribution')
    ax.set_title(image_name, fontsize='x-large')
    plt.legend(loc='best')
    plot_out_dir = create_dir(out_dir, 'distribution_plots')
    out_name = os.path.join(plot_out_dir, '%s.png' % image_name)
    plt.savefig(out_name)
    plt.close(fig)

def train_step(real_sample):
    for i in range(8):
        with tf.GradientTape() as disc_tape:
            fake_sample, dist_cdf = generator()
            fake_output = discriminator(fake_sample)
            real_output = discriminator(real_sample)
            disc_loss = discriminator_loss(real_output, fake_output)
        disc_gradients = disc_tape.gradient(disc_loss, discriminator.trainable_variables)
        discriminator_optimizer.apply_gradients(zip(disc_gradients, discriminator.trainable_variables))
    with tf.GradientTape() as gen_tape:
        fake_sample, dist_cdf = generator()
        fake_output = discriminator(fake_sample)
        dist_loss, gen_loss = generator_loss(fake_output, dist_cdf)
        total_loss = 0.995 * dist_loss + 0.005 * gen_loss
    gen_gradients = gen_tape.gradient(total_loss, [generator.lamb, generator.mu, generator.sigma])
    generator_optimizer.apply_gradients(zip(gen_gradients, [generator.lamb, generator.mu, generator.sigma]))
    return disc_loss, dist_loss, gen_loss

log_dir = 'logs'
summary_writer = tf.summary.create_file_writer(log_dir)

def create_dir(parent_dir, dir_name):
    """
    Create directory.
    :param parent_dir: parent directory
    :param dir_name: new directory name
    :return: new directory path
    """
    out_dir = os.path.join(parent_dir, dir_name)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    return str(os.path.join(parent_dir, dir_name))

def train(FACS_data):
    num_iter = 0
    min_dist_loss = 10
    num_iter_list = []
    dist_loss_list = []
    gen_loss_list = []
    disc_loss_list = []
    for real_sample in FACS_data:
        num_iter += 1
        disc_loss, dist_loss, gen_loss = train_step(real_sample)
        num_iter_list.append(num_iter)
        dist_loss_list.append(dist_loss)
        gen_loss_list.append(gen_loss)
        disc_loss_list.append(disc_loss)
        mu = generator.mu.numpy()
        sigma = generator.sigma.numpy()
        lamb = generator.lamb.numpy()
        if num_iter % 50 == 0:
            with summary_writer.as_default():
                tf.summary.scalar('discriminator_loss', disc_loss, step=num_iter)
                tf.summary.scalar('distribution_loss', dist_loss, step=num_iter)
                tf.summary.scalar('generator_loss', gen_loss, step=num_iter)
            # Plot training curve
            fig, ax = plt.subplots()
            ax.plot(num_iter_list, dist_loss_list, label='distribution_loss')
            ax.plot(num_iter_list, gen_loss_list, label='generator_loss')
            ax.plot(num_iter_list, disc_loss_list, label='discriminator_loss')
            ax.legend()
            ax.set_xlabel('num_iter')
            ax.set_ylabel('loss')
            ax.set_title('Training curve')
            training_out_dir = create_dir(out_dir, 'training_curve')
            out_name = os.path.join(training_out_dir, "num_iter%d.png" % num_iter)
            plt.savefig(out_name)
            plt.close(fig)

        print('num_iter:%d, disc_loss:%f, dist_loss:%f, gen_loss:%f' % (num_iter, disc_loss, dist_loss, gen_loss))

        distribution_plot(lamb, mu, sigma, image_name='num_iter%d' % num_iter)
        if dist_loss < min_dist_loss:
            min_dist_loss = dist_loss
            with open(os.path.join(out_dir, 'mu.pickle'), 'wb') as f:
                pickle.dump(mu, f)
            with open(os.path.join(out_dir, 'sigma.pickle'), 'wb') as f:
                pickle.dump(sigma, f)
            with open(os.path.join(out_dir, 'lamb.pickle'), 'wb') as f:
                pickle.dump(lamb, f)

train(FACS_data)

tf.keras.backend.clear_session()
gc.collect()