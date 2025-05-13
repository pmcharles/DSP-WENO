# Created by: Deep Ray, University of Maryland
# Date: March 25, 2023
# Modified by: Phil Charles (8/24/23) for strong sign property constraint

import argparse, textwrap

formatter = lambda prog: argparse.HelpFormatter(prog,max_help_position=50)

def cla():
    parser = argparse.ArgumentParser(description='list of arguments',formatter_class=formatter)

    # Data parameters
    parser.add_argument('--data_file' , type=str, default='DSP-WENO_training_data.dat', help=textwrap.dedent('''Data file containing cell values and interface values'''))
    parser.add_argument('--tanh_theta', type=bool, default=True, help=textwrap.dedent('''Whether or not we are feeding in tanh(theta)'''))
    parser.add_argument('--n_train'   , type=float, default=0.6, help=textwrap.dedent('''Fraction of data to be used as training samples.'''))
    parser.add_argument('--n_val'     , type=float, default=0.2, help=textwrap.dedent('''Fraction of data to be used as validation samples.'''))
    parser.add_argument('--n_test'    , type=float, default=0.2, help=textwrap.dedent('''Fraction of data to be used as test samples.'''))
    parser.add_argument('--scaling'   , type=bool, default=False, help=textwrap.dedent('''Whether or not to scale the data.'''))
    parser.add_argument('--scale_met' , type=int, default=1, help=textwrap.dedent('''Scaling method used.'''))

    # Network parameters
    parser.add_argument('--h_widths'     , nargs='+', type=int, default=[5,5,5], help=textwrap.dedent('''Width of hidden layers of network'''))
    parser.add_argument('--act_func'     , type=str, default='ReLU', help=textwrap.dedent('''Activation function'''))
    parser.add_argument('--act_func_par' , type=float, default=0.2, help=textwrap.dedent('''Only applicable if using LeakyReLU'''))
    parser.add_argument('--lr'           , type=float, default=1e-3, help=textwrap.dedent('''Learning rate'''))
    parser.add_argument('--max_epoch'    , type=int, default=50, help=textwrap.dedent('''Maximum number of epochs'''))
    parser.add_argument('--batch_size'   , type=int, default=500, help=textwrap.dedent('''Batch size while training'''))
    parser.add_argument('--loss_func'    , type=str, default='MSE', help=textwrap.dedent('''Loss function'''))
    parser.add_argument('--adam_beta1'   , type=float, default=0.5, help=textwrap.dedent('''Beta 1 coefficient for Adam optimizer'''))
    parser.add_argument('--adam_beta2'   , type=float, default=0.9, help=textwrap.dedent('''Beta 2 coefficient for Adam optimizer'''))
    parser.add_argument('--weight_dec'   , type=float, default=1e-5, help=textwrap.dedent('''Weight decay regularizer'''))
    parser.add_argument('--seed_no'      , type=int, default=6008, help=textwrap.dedent('''Set the random seed'''))
    parser.add_argument('--input_size'   , type=int, default=15, help=textwrap.dedent('''How many input values per data point'''))
    
    # Output parameters
    parser.add_argument('--val_freq'     , type=int, default=2, help=textwrap.dedent('''Number of epochs after validation loss is computed'''))
    parser.add_argument('--base_dir'     , type=str, default='DSP-WENO_trained_models', help=textwrap.dedent('''Base directory where all trained models and results are saved'''))

    
    return parser.parse_args()