import os
import sys
from pathlib import Path
import argparse
import numpy as np
import nibabel as nib

import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision.models as models
import torch.optim as optim
from torch.utils import data

## Please specify the directory of the LPCNN code from Github
# LPCNN_HOME=Path('/project/3015069.05/bids/code/LPCNN/')

sys.path.append(os.path.join(LPCNN_HOME,'LPCNN'))

from lib.utils import *
from lib.tool.tool import qsm_display

os.environ['CUDA_VISIBLE_DEVICES'] = '0'

# output directory
root_dir = os.path.join(LPCNN_HOME,'data')
vis_output_path = Path(os.path.join('LPCNN','test_result'))

# data statistics
gamma = 42.57747892
gt_mean = np.load(os.path.join(root_dir,'numpy_data','whole','list','train_gt_mean.npy'))
gt_std = np.load(os.path.join(root_dir,'numpy_data','whole','list','train_gt_std.npy'))


def main(args):
	
	device = torch.device('cuda:0' if not args.no_cuda else 'cpu')

	# name
	example_name = args.model_arch + args.save_name
	
	# load data
	phase_path_list = []
	with open(str(Path(args.phase_file))) as f:
		for line in f:
			phase_path_list.append(line.strip('\n'))
	x, y, z = nib.load(str(Path(phase_path_list[0]))).get_fdata().shape
	phase_data_list = np.zeros((x, y, z, args.number))

	phase_example = str(Path(phase_path_list[0]))
	for i in range(args.number):
		phase_data_list[:, :, :, i] = nib.load(str(Path(phase_path_list[i]))).get_fdata() / (args.tesla * gamma)

	dipole_data_list = ''
	with open(str(Path(args.dipole_file))) as f:
		for line in f:
			dipole_data_list = dipole_data_list + str(Path(line.strip('\n'))) + ' '
	
	with open(str(Path(args.mask_file))) as f:
		for line in f:
			mask_data = nib.load(str(Path(line.strip('\n')))).get_fdata()

	if not args.gt_file == None:
		with open(str(Path(args.gt_file))) as f:
			for line in f:
				gt_data = nib.load(str(Path(line.strip('\n')))).get_fdata()
	
	# crop
	crop = args.crop
	if crop == None:
		crop = [0, 0, 0, 0, 0, 0]
	if len(crop) != 6:
		print('ERROR: crop parameters should be 6 intgers')
		exit()

	# load model
	model = chooseModel(args, root_dir / Path(os.path.join('numpy_data','whole','list')))
	model.load_state_dict(torch.load(args.resume_file, map_location=device)['model_state'])
	model.to(device)
	print(args.model_arch + ' loaded.')
	
	# parallel model
	if args.gpu_num > 1:
		model = nn.DataParallel(model)

	model_name = Path(args.resume_file).parts[-1].split('.')[0]
	print(model_name)

	mse_loss = 0
	ssim_perf = 0
	psnr_perf = 0

	model.eval()

	with torch.no_grad():
	
		#cuda
		phase_data_list = torch.from_numpy(phase_data_list[np.newaxis,np.newaxis,crop[0]:x-crop[1],crop[2]:y-crop[3],crop[4]:z-crop[5],:])
		phase_data_list = phase_data_list.to(device, dtype=torch.float)

		mask_data = torch.from_numpy(mask_data[np.newaxis,crop[0]:x-crop[1],crop[2]:y-crop[3],crop[4]:z-crop[5]])
		mask_data = mask_data.to(device, dtype=torch.float)
	
		output_data = model(phase_data_list, (dipole_data_list,), mask_data.unsqueeze(1))

		mask = torch.squeeze(mask_data, 0).cpu().numpy()[:,:,:,np.newaxis]

		#denormalize
		if args.norm:

			og_output = torch.squeeze(output_data, 0).permute(1, 2, 3, 0).cpu().numpy() * gt_std
			og_output = og_output + gt_mean
			og_output = og_output * mask

		
		if not args.gt_file == None:
			og_gt = gt_data[:, :, :, np.newaxis] * mask
			mse_loss += np.sqrt(qsm_mse(og_gt, og_output, mask, roi=True))
			ssim_perf += qsm_ssim(og_gt, og_output, mask, root_dir)	
			psnr_perf += qsm_psnr(og_gt, og_output, mask, root_dir, roi=True)

		output_path = vis_output_path / model_name
		output_path.mkdir(parents=True, exist_ok=True)

		save_name = output_path / example_name
	
		og_output = np.pad(og_output, ((crop[0],crop[1]), (crop[2],crop[3]), (crop[4],crop[5]), (0,0)), mode='constant')

		qsm_display(og_output, phase_example, torch.squeeze(mask_data, 0).cpu().numpy(), out_name=str(save_name))

	if not args.gt_file == None:
		avg_mse_loss = mse_loss / np.sqrt(whole_validation_mse)
		avg_ssim_perf = ssim_perf
		avg_psnr_perf = psnr_perf

		print('##Test RMSE: %.8f PSNR: %.8f SSIM: %.8f' %(avg_mse_loss, avg_psnr_perf, avg_ssim_perf))


parser = argparse.ArgumentParser(description='QSM Inference')
parser.add_argument('--save_name', type=str, default='_example', help='output name')
parser.add_argument('--number', type=int, default=3, choices=[1, 2, 3], help='input phase number')
parser.add_argument('--phase_file', type=str, default='test_data/one/phase_data1.txt', help='the testing phase data')
parser.add_argument('--dipole_file', type=str, default='test_data/one/dipole_data1.txt', help='the testing dipole data')
parser.add_argument('--mask_file', type=str, default='test_data/one/mask_data1.txt', help='the testing mask data')
parser.add_argument('--gt_file', type=str, default=None, help='the testing gt data')
parser.add_argument('--tesla', default=7, type=int, choices=[3, 7], help='B0 tesla(default: 7)')
parser.add_argument('--norm', action='store_true', default=True, help='normalize data')
parser.add_argument('--gpu_num', default=1, type=int, choices=[1, 2, 3, 4], help='number of gpu (default: 1)')
parser.add_argument('--model_arch', default='lpcnn', choices=['lpcnn'], help='network model (default: lpcnn)')
parser.add_argument('--no_cuda', action='store_true', default=False, help='disables CUDA training')
parser.add_argument('--resume_file', type=str, default=None, help='the checkpoint file to resume from')

parser.add_argument('--crop', type=int, nargs='+',  help='crops redundant margin')

if __name__ == '__main__':
	args = parser.parse_args()
	main(args)
