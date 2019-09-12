# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with JabberDock;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# Please reference this software and its authors in any works that make use of it
# Author : Lucas Rudden, lucas.rudden1@btinternet.com

import numpy as np

def create_pqr_dat(in_file, out_file, template="../biobox/classes/amber03.dat"):
    """
    Convert an rtp file (e.g. in gromacs top folder) into an appropiate pqr format that can be interpreted by biobox to form the
    STID maps

    :param in_file: Name of input rtp file, e.g. aminoacids.rtp in amber03.ff/ folder 
    :param out_file: Name of output file that can be read by the pqr interpreter (import_pqr in biobox.Molecule)
    :param template: Template of pqr file. We recommend just using the amber03.dat shipped with biobox
    """

    data = []
    atoms = False

    input_file = open(in_file, "r")
    prevline1 = ""
    prevline2 = "" # store the previous line so we have the residue too

    for line in input_file:
        if line.strip() == "[ atoms ]":
            atoms = True # now we're in the atoms regime (always occurs before atoms are listed)
            prevline2 = prevline1
        elif atoms and len(line.split()) == 4:
            tmp = np.insert(arr=line.split(), obj=0, values=prevline2.split(' ')[1])
            data.append([tmp[0], tmp[1], tmp[3], '', tmp[2]])
        else:
            atoms = False
        prevline1 = line

    # now we want to import the data file that we're using as a template
    # This gives us the format for the pqr filetype. You can use older .dat files (e.g. amber14sb) for this
    tem_file = np.loadtxt(template, dtype=str)
    uni = np.unique(cp_file[:, 4])
    data = np.asarray(data)

    for i in uni:
        vdw_place = cp_file[:, 4] == i
        # all the vdw values are the same anyway for the atomtype
        vdw_val = cp_file[vdw_place][:, 3][0]
        atomtypes = data[:,4] == i
        data[atomtypes,3] = vdw_val
    
    f_out = open(out_file, "w")
    for i in range(0, np.shape(data)[0]):
        L = "%-3s     %-8s%-16s%-8s%s\n"%(str(data[i,0]), str(data[i,1]), str(data[i,2]), str(data[i,3]), str(data[i,4]))
        f_out.write(L)
    f_out.close()
 
