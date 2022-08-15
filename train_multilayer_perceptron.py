# ----------------------------------------------------------------------------
# train_multilayer_perceptron.py 
# Author: Antoine Beauchamp
# Created: December 15th, 2020

"""
Train a multi-layer perceptron

Description
-----------

"""

# Packages -------------------------------------------------------------------

import argparse
import os
import random
import csv
import pandas                 as pd
import numpy                  as np
from datatable                import fread

import torch
import torch.nn.functional    as F
from torch                    import nn
from torch.optim              import SGD, AdamW
from torch.optim.lr_scheduler import OneCycleLR
from torch.cuda               import is_available

from skorch                   import NeuralNetClassifier
from skorch.callbacks         import LRScheduler
from skorch.helper            import DataFrameTransformer

from sklearn.metrics          import accuracy_score, confusion_matrix


# Functions ------------------------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        "--outdir",
        type = str,
        default = 'data/MLP_outcomes/',
        help = "Directory in which to write neural net outcomes."
    )
    
    parser.add_argument(
        '--training',
        type = str,
        help = "Path to CSV file containing training data."
    )
    
    parser.add_argument(
        '--mousetransform',
        type = str,
        help = ("Path to CSV file containing mouse data to transform")
    )
    
    parser.add_argument(
        '--humantransform',
        type = str,
        help = ("Path to CSV file containing human data to transform")
    )
    
    parser.add_argument(
        "--nunits",
        type = int,
        default = 500,
        help = "Number of hidden units in the network."
    )
    
    parser.add_argument(
        "--L2",
        type = float,
        default = 1e-6,
        help = "Weight decay"
    )
    
    parser.add_argument(
        "--nepochs",
        type = int,
        default = 200,
        help = "Number of epochs to train over."
    )
    
    parser.add_argument(
        '--totalsteps',
        type = int,
        help = "Number of steps to use in optimizer."
    )
    
    parser.add_argument(
        "--learningrate",
        type = float,
        default = 1e-5,
        help = "Learning rate during training."
    )
    
     parser.add_argument(
        '--optimizer',
        type = str,
        default = 'SGD'
    )
    
    parser.add_argument(
        "--confusionmatrix",
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = ("Option to compute confusion matrix from training set "
                "predictions.")
    )
    
    parser.add_argument(
        '--seed',
        type = int,
        help = ("Random seed")
    )
    
    parser.add_argument(
        '--saveparams',
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = ("Save neural network parameters?")
    )
    
    parser.add_argument(
        '--paramsheader',
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = ("Save neural network parameters header?")
    )
    
    parser.add_argument(
        '--verbose',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = ("Verbosity.")
    )
        
    args = vars(parser.parse_args())
    
    return args


def prepare_data(data, labelcol = None):

    """
    """

    if labelcol is not None:
    
        ind_labels = data.columns.str.match(labelcol)

        input_data = data.loc[:, ~ind_labels].copy()
        label_data = data.loc[:, ind_labels].copy()

        label_data.loc[:, labelcol] = label_data.loc[:, labelcol].astype('category')

        y = DataFrameTransformer().fit_transform(label_data)[labelcol]
        
    else:
        
        input_data = data.copy()
        y = None

    X = DataFrameTransformer().fit_transform(input_data)['X']
        
    return X, y


def define_classifier(input_units, output_units, hidden_units, weight_decay, max_epochs, total_steps, learning_rate, optimizer, device):

    """
    """
    
    #Define network architecture
    class ClassifierModule(nn.Module):
        def __init__(
            self,
            input_units,
            output_units,
            hidden_units,
            apply_output_layer = True #Flag to apply output layer
        ):
            super(ClassifierModule, self).__init__()

            self.apply_output_layer = apply_output_layer

            self.hidden1 = nn.Linear(input_units, hidden_units)
            self.hidden2 = nn.Linear(hidden_units, hidden_units)
            self.hidden3 = nn.Linear(hidden_units, hidden_units)
            self.output = nn.Linear(hidden_units, output_units)

        def forward(self, X, **kwargs):
            X = F.relu(self.hidden1(X))
            X = F.relu(self.hidden2(X))
            X = F.relu(self.hidden3(X))

            #If flag is True, apply output layer
            if self.apply_output_layer is True:
                X = F.softmax(self.output(X), dim = -1)

            return X

    #Create the classifier
    net = NeuralNetClassifier(
        ClassifierModule(input_units = input_units,
                         output_units = output_units,
                         hidden_units = hidden_units),
        train_split = None,
        optimizer = optimizer,
        optimizer__weight_decay = weight_decay,
        max_epochs = max_epochs,
        callbacks = [('lr_scheduler',
                      LRScheduler(policy=OneCycleLR,
                                  total_steps=total_steps,
                                  cycle_momentum=False,  
                                  max_lr=learning_rate))],
        device = device)

    return net


def calc_confusion_matrix(data, y_true, y_pred, labelcol):

    """
    """

    df_confmat = pd.DataFrame(confusion_matrix(y_true, y_pred))

    df_labels = data[[labelcol]].copy()
    df_labels['Dummy'] = y_true
    df_labels = df_labels.sort_values('Dummy').drop_duplicates().reset_index(drop = True)

    df_confmat.columns = df_labels[labelcol].astype('str')

    df_confmat['TrueLabels'] = (df_labels[labelcol]
                                .astype('str')
                                .reset_index(drop = True))

    return df_confmat
    
    
def transform_input_space(data, network, labelcol = None):

    """
    """
    
    X, y = prepare_data(data = data, 
                        labelcol = labelcol)

    df_transformed = pd.DataFrame(network.predict_proba(X))

    if labelcol is not None:
        df_transformed[labelcol] = data[labelcol]

    return df_transformed
    
    
# Main ------------------------------------------------------------------------

def main():

    #Load command line arguments
    args = parse_args()
    
    outdir = args['outdir']
    confmat = True if args['confusionmatrix'] == 'true' else False
    saveparams = True if args['saveparams'] == 'true' else False
    paramsheader = True if args['paramsheader'] == 'true' else False
    verbose = True if args['verbose'] == 'true' else False
    
    outdir = os.path.join(outdir, '')
    if os.path.exists(outdir) == False:
        print('Output directory {} not found. Creating it...'.format(outdir))
        os.mkdir(outdir)
    
    
    # Import data -------------------------------------------------------------

    if verbose:
        print("Importing training data...")
    
    training_file = args['training']
    df_training = fread(training_file, header = True).to_pandas()

    
    # Process data ------------------------------------------------------------

    if verbose:
        print("Preparing data for learning...")
    
    X, y = prepare_data(data = df_training, 
                        labelcol = 'Region')
    

    # Train the network ------------------------------------------------------

    if verbose:
        print("Initializing neural network...")

    seed = args['seed']
    if seed is not None:
        np.random.seed(seed)
        torch.manual_seed(seed)
        random.seed(seed)
        
    #Get network parameters from command line args
    hidden_units = args['nunits']
    weight_decay = args['L2']
    max_epochs = args['nepochs']
    total_steps = args['totalsteps']
    learning_rate = args['learningrate']
    optimizer = args['optimizer']
    nlabels = len(np.unique(y))
    
    if total_steps is None:
        total_steps = max_epochs
        
    if optimizer == 'AdamW':
        optimizer = AdamW
    elif optimizer == 'SGD':
        optimizer = SGD
    else:
        raise ValueError
    
    if is_available() == True:
        print("GPU available. Training network using GPU...")
        device = 'cuda'
    else:
        print("GPU unavailable. Training network using CPU...")
        device = 'cpu'
    
    #Define the model object
    net = define_classifier(input_units = X.shape[1],
                            output_units = nlabels,
                            hidden_units = hidden_units,
                            weight_decay = weight_decay,
                            max_epochs = max_epochs,
                            total_steps = total_steps,
                            learning_rate = learning_rate,
                            optimizer = optimizer,
                            device = device)
    
    #Fit the network
    net.fit(X, y)
    
    #Predict training labels
    y_pred = net.predict(X)
    
    #Save parameters
    if saveparams:
    
        loss_func = net.__dict__['criterion'].__name__
        train_loss = net.__dict__['history_'][max_epochs-1]['train_loss']
        train_accuracy = accuracy_score(y, y_pred)
        dur_list = [net.__dict__['history_'][i]['dur'] 
                    for i in range(len(net.__dict__['history_']))]
        dur_mean = np.mean(dur_list)
    
        params = [hidden_units,
                  weight_decay,
                  max_epochs,
                  total_steps,
                  learning_rate,
                  seed,
                  nlabels,
                  loss_func, 
                  train_loss, 
                  train_accuracy,
                  dur_mean]
        
        if paramsheader:
            
            header = ['hidden_units',
                      'weight_decay',
                      'nepochs',
                      'total_steps',
                      'learning_rate',
                      'random_seed',
                      'nlabels',
                      'loss_function',
                      'training_loss',
                      'training_accuracy',
                      'mean_epoch_duration']
            
        file_params = os.path.join(outdir, 'MLP_params.csv')

        with open(file_params, 'w') as file:
            
            write = csv.writer(file)
            
            if paramsheader:
                write.writerow(header)
                
            write.writerow(params)
    
    
    # Compute training confusion matrix --------------------------------------
 
    if confmat:
        
        if verbose:
            print("Computing confusion matrix from training set...")
        
        df_confmat = calc_confusion_matrix(data = df_training,
                                           y_true = y,
                                           y_pred = y_pred,
                                           labelcol = 'Region')
        
               
        if verbose:
            print("Writing the confusion matrix to file...")
        
        #File to save confusion matrix
        file_confmat = "MLP"+\
        "_labels"+str(nlabels)+\
        "_layers3"+\
        "_units"+str(args['nunits'])+\
        "_L2"+str(args['L2'])+\
        "_confusionmatrix_training.csv"
        
        #Write confusion matrix to file
        df_confmat.to_csv(os.path.join(outdir, file_confmat), 
                          index = False)


    # Predict label probabilities for mouse/human data -----------------------

    mouse_transform = args['mousetransform']
    human_transform = args['humantransform']

    #If no mouse file given, use training data
    if mouse_transform is not None:
        
        if verbose:
            print("Importing mouse data to transform...")
        
        df_mouse = fread(mouse_transform, header = True).to_pandas()
        
    else: 
        df_mouse = df_training
        
    if 'Region' in df_mouse.columns:
        mouselabels = 'Region'
    else:
        mouselabels = None

    if verbose:
        print("Computing mouse label probabilities...")
        
    #Compute mouse probabilities
    df_mouse_prob = transform_input_space(data = df_mouse,
                                          network = net,
                                          labelcol = mouselabels)
    
    #Assign column names 
    if mouselabels is None:
        df_mouse_prob.columns = np.unique(df_training['Region'])
    else:
        df_mouse_prob.columns = np.append(np.unique(df_training['Region']), mouselabels)
        
    if verbose:
        print("Writing mouse probabilities to file...")

    #File to save mouse probabilities
    file_mouse_prob = 'MLP'+\
    '_labels'+str(nlabels)+\
    '_layers3'+\
    '_units'+str(args['nunits'])+\
    '_L2'+str(args['L2'])+\
    '_mouseprob'
    
    if mouselabels is not None:
        #Number of labels in the mouse file
        nlabels_mouse = len(np.unique(df_mouse[mouselabels]))
        file_mouse_prob = file_mouse_prob+'_'+str(nlabels_mouse)
        
    file_mouse_prob = file_mouse_prob+'.csv'

    #Write probability data to file
    df_mouse_prob.to_csv(os.path.join(outdir, file_mouse_prob),
                         index = False)
    
    if human_transform is not None:
        
        if verbose:
            print("Importing human data to transform...")
        
        df_human = (fread(human_transform, header = True)
                    .to_pandas())
        
        if 'Region' in df_human.columns:
            humanlabels = 'Region'
        else:
            humanlabels = None
        
        if verbose:
            print("Computing human label probabilities...")
        
        df_human_prob = transform_input_space(data = df_human,
                                              network = net,
                                              labelcol = humanlabels)
        
        #Assign column names 
        if humanlabels is None:
            df_human_prob.columns = np.unique(df_training['Region'])
        else:
            df_human_prob.columns = np.append(np.unique(df_training['Region']), humanlabels)
        
        if verbose:
            print("Writing human probabilities to file...")
        
        file_human_prob = 'MLP'+\
           '_labels'+str(nlabels)+\
        '_layers3'+\
        '_units'+str(args['nunits'])+\
        '_L2'+str(args['L2'])+\
        '_humanprob'
        
        if humanlabels is not None:
            #Number of labels in the human file
            nlabels_human = len(np.unique(df_human[humanlabels]))
            file_human_prob = file_human_prob+'_'+str(nlabels_human)
            
        file_human_prob = file_human_prob+'.csv'
        
        df_human_prob.to_csv(os.path.join(outdir, file_human_prob),
                             index = False)
    
       
    # Transform the mouse/human data to the latent space ------------------------------

    #Change the mode of the network to extract a hidden layer
    net.module_.apply_output_layer = False
    
    if verbose:
        print("Transforming mouse data to latent space...")
    
    df_mouse_transform = transform_input_space(data = df_mouse,
                                               network = net,
                                               labelcol = mouselabels)
    
    if verbose:
        print("Writing mouse latent space data to file...")
    
    file_mouse_transform = 'MLP'+\
    '_labels'+str(nlabels)+\
    '_layers3'+\
    '_units'+str(args['nunits'])+\
    '_L2'+str(args['L2'])+\
    '_mousetransform'

    if mouselabels is not None:
        file_mouse_transform = file_mouse_transform+'_'+str(nlabels_mouse)
        
    file_mouse_transform = file_mouse_transform+'.csv'
    
    df_mouse_transform.to_csv(os.path.join(outdir, file_mouse_transform), 
                              index = False)
    
    if human_transform is not None:
        
        if verbose:
            print("Transforming human data to latent space...")
        
        df_human_transform = transform_input_space(data = df_human,
                                                   network = net,
                                                   labelcol = humanlabels)
        
        if verbose:
            print("Writing human latent space data to file...")
        
        file_human_transform = 'MLP'+\
        '_labels'+str(nlabels)+\
        '_layers3'+\
        '_units'+str(args['nunits'])+\
        '_L2'+str(args['L2'])+\
        '_humantransform'
        
        if humanlabels is not None:
            file_human_transform = file_human_transform+'_'+str(nlabels_human)
            
        file_human_transform = file_human_transform+'.csv'
        
        df_human_transform.to_csv(os.path.join(outdir, file_human_transform),
                                  index = False)
    
        
if __name__ == "__main__":
    main()
