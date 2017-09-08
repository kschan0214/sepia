The FAst Nonlinear Susceptibility Inversion (FANSI) Toolbox was
developed by Carlos Milovic at the Biomedical Imaging Center at 
Pontificia Universidad Católica de Chile, in 2017, with the 
colaboration of:
Berkin Bilgic and Bo Zhao at Martinos Center for Biomedical 
Imaging, Harvard Medical School, MA, USA
Julio Acosta-Cabronero at Wellcome Trust Centre for Neuroimaging, 
Institute of Neurology, University College London, London, UK, 
and German Center for Neurodegenerative Diseases (DZNE), 
Magdeburg, Germany
and Cristian Tejos at Department of Electrical Engineering, 
Pontificia Universidad Catolica de Chile, Santiago, Chile and 
the Biomedical Imaging Center at Pontificia Universidad Católica 
de Chile.
 
Please cite the following publications if you use this toolbox or any 
code in it:


Bilgic B, Fan AP, Polimeni JR, Cauley SF, Bianciardi M, 
Adalsteinsson E, Wald LL, Setsompop K. Fast quantitative susceptibility 
mapping with L1-regularization and automatic parameter selection. 
Magn Reson Med. 2014;72:1444-1459.

Bilgic B., Chatnuntawech I., Langkammer C., Setsompop K.; Sparse Methods 
for Quantitative Susceptibility Mapping; Wavelets and Sparsity XVI, 
SPIE 2015

Milovic C, Bilgic B, Zhao B, Acosta-Cabronero J, Tejos C. A Fast 
Algorithm for Nonlinear QSM Reconstruction with Variational Penalties. 
Proceedings of the Fourth International Workshop on MRI Phase Contrast 
& Quantitative Susceptibility Mapping. 2016;1:132.

Milovic C, Bilgic B, Zhao B, Acosta-Cabronero J, Tejos C. Fast Nonlinear 
QSM Method with Variational Regularizations. Proceedings of the 
ISMRM 2017.
[UPDATE WHEN AVALAIBLE WITH PAPER]

 
 
 
 
The code of this toolbox was based on the source code released by 
Berkin Bilgic at http://martinos.org/~berkin/software.html


We provide an analytic Susceptibility brain phantom based on the 
one developed by Christian Langkammer, et al in "Fast quantitative 
susceptibility mapping using 3D EPI and total generalized variation",
and C. Wisnieff, et al in "Magnetic susceptibility anisotropy: 
cylindrical symmetry from macroscopically ordered anisotropic 
molecules and accuracy of MRI measurements using few orientations." 
Please cite our work and theyrs if you use this phantom.
 
The functions used to calculate quality metrics are based on the source 
code provided at the 4th International Workshop on MRI Phase Contrast 
& Quantitative Susceptibility Mapping, September 26th - 28th 2016, at
Medical University of Graz, Austria, for its QSM Challenge. 
Please cite the Challenge report paper if you use this code and metrics.
[UPDATE WHEN AVALAIBLE]
 
The function to calculate the SSIM quality index uses code 
by Zhou Wang. Please see the header information for more details.
Modified by Berkin Bilgic for QSM data.
If you find this useful, please cite:
Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
quality assessment: From error measurement to structural similarity"
IEEE Transactios on Image Processing, vol. 13, no. 4, Apr. 2004.
 
The function to calculate the Mutual Information quality index uses code 
by R. Moddemeijer, at http://www.cs.rug.nl/~rudy/matlab/
If you find this useful, please cite:
Moddemeijer, R. On Estimation of Entropy and Mutual Information of 
Continuous Distributions, Signal Processing, 1989, vol. 16, nr. 3, 
pp. 233-246
 
 
 
DISCLAIMER:
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
SHALL CARLOS MILOVIC OR HIS CONTRIBUTORS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY 
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
BUSINESS INTERRUPTION; PROCUREMENT OF SUBSTITUTE GOODS 
OR SERVICES; AND LOSS OF USE, DATA OR PROFITS) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
DAMAGE.

