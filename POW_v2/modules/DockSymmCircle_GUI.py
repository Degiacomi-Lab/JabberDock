#! /usr/bin/env python
# Copyright (c) 2012 EPFL (Ecole Polytechnique federale de Lausanne)
# Laboratory for Biomolecular Modeling, School of Life Sciences
#
# POW is free software ;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation ;
# either version 2 of the License, or (at your option) any later version.
# POW is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with POW ;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# Author : Matteo Degiacomi, matteothomas.degiacomi@epfl.ch
# Web site : http://lbm.epfl.ch


import os, sys
import wx
import subprocess
import DockSymmCircle as module

class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        self.dirname=''
        # A "-1" in the size parameter instructs wxWidgets to use the default size.
        # In this case, we select 200px width and the default height.
        wx.Frame.__init__(self, parent, title=title, size=(600,500))
        self.CreateStatusBar() # A Status bar in the bottom of the window


        self.SetBackgroundColour("white")

        # Setting up the menu bar
        filemenu= wx.Menu()
        menuLoad = filemenu.Append(wx.ID_OPEN, "&Load..."," Load a setup file")
        menuSave = filemenu.Append(wx.ID_SAVEAS, "&Save as..."," Save the current setup")
        menuHelp = filemenu.Append(wx.ID_HELP, "&Help"," Optimization user manual")
        menuAbout= filemenu.Append(wx.ID_ABOUT, "&About"," Information about this program")
        filemenu.AppendSeparator()
        menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")

        # Creating the menu bar
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

        # launch and verify buttons
        self.launch=wx.Button(self, label="launch!")
        self.verify=wx.Button(self, label="verify")
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        buttons.Add(self.verify, 0 )
        buttons.Add(self.launch, 0 )

        #panels
        self.I=IOPanel(self)
        self.C=ConstraintPanel(self)
        self.B=BoundaryPanel(self)
        self.R=ReceptorPanel(self)
        self.S=ProtocolPanel(self)

        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.I, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.B, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.R, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.C, 1, wx.EXPAND | wx.ALL, 3 )
        #box.Add(self.P, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.S, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(buttons, 0, wx.ALIGN_CENTER)
        self.SetSizer(box)
        self.SetSize((400,800))
        self.Centre()

        # bind events to menu items
        self.Bind(wx.EVT_MENU, self.OnLoad, menuLoad)
        self.Bind(wx.EVT_MENU, self.OnSave, menuSave)
        self.Bind(wx.EVT_MENU, self.OnHelp, menuHelp)
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
        self.Bind(wx.EVT_BUTTON, self.OnLaunch,self.launch)
        self.Bind(wx.EVT_BUTTON, self.OnVerify,self.verify)

        self.Show()


    def OnAbout(self,e):
        # Create a message dialog box
        dlg = wx.MessageDialog(self, " Prediction of circular protein assemblies@POW \n Matteo Degiacomi \n 2012", "About POW", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.

    def OnExit(self,e):
        """ Close the interface"""
        self.Close(True)  # Close the frame.

    def OnSave(self,e):
        """ Save a file"""
        dlg = wx.FileDialog(self, "Chose a name for your setup file!", self.dirname, "", "*.*", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.makeFile(os.path.join(self.dirname, self.filename))
        dlg.Destroy()

    def OnHelp(self,e):
        """ Open lbm webpage """
        dlg = wx.MessageDialog(self, "Code updates and user manual available at \n lbm.epf.ch/resources", "POW help", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.


    def OnLoad(self,e):
        """ Open a file"""
        infile=-1
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            infile =os.path.join(self.dirname, self.filename)
        dlg.Destroy()

        #call file parsing
        if infile!=-1:
            if os.path.isfile(infile):
                self.parse(infile)

    def init(self):
        """ Reset all the form fields"""
        self.S.editsteps.SetValue('')
        self.S.editpart.SetValue('')
        self.S.editrep.SetValue('')
        self.C.setconstraint.SetValue('')
        for __ in xrange(0,self.C.lc.GetItemCount(),1):
            self.C.lc.DeleteItem()

    def parse(self,infile):
        """ Input file parser"""

        self.init()

        params=module.Parser() #read user defined variables
        params.add_standard() #add default variables
        params.set_default_values() #set default values to all defined variables
        params.parse(infile) #parse input file

        #protocol panel
        self.S.editrep.SetValue(str(params.repeat))
        self.S.editsteps.SetValue(str(params.max_steps))
        self.S.editpart.SetValue(str(params.n_particles))

        #protocol advanced settings
        self.S.editneighsize=str(params.neigh_size)
        self.S.editneigh=str(params.neigh_type)
        self.S.editinertiamax=float(params.inertia_max)
        self.S.editinertiamin=float(params.inertia_min)
        self.S.editpersonal=float(params.cp)
        self.S.editglobal=float(params.cn)

        #constraint panel
        #insert data from target values list
        for x in xrange(0,len(params.target),1):
            self.C.lc.InsertStringItem(x, str(x+1))
            if str(params.target[x]) != "NaN":
                self.C.lc.SetStringItem(x, 1, str(params.target[x]))
            #if len(params.keyword)!=0 and str(params.keyword(x)) != "NaN":
            #    self.B.lc.SetStringItem(x, 2, str(params.keyword(x)))

        self.C.setconstraint.SetValue(str(params.constraint))
        self.C.editmix.SetValue(str(params.mix_weight))
        self.C.editfilter.SetValue(str(params.accept))
        #self.add('detectClash','detect_clash','str',"on")

        #input panel
        self.I.setdegree.SetValue(str(params.degree))
        self.I.editassemblystyle.SetValue(str(params.style))
        self.I.setstructure.SetValue(str(params.pdb_file_name))
        self.I.settop.SetValue(str(params.topology))
        self.I.setproj.SetValue(str(params.proj_file))
        self.I.setratio.SetValue(str(params.ratio))
        #self.add('align','align','str',"no")
        #self.add('trajectory','trajectory','str',"NA")
        #self.add('trajSelection','trajselection','str',"NA")

        print str(params.pdb_file_name)

        #boundary panel
        if len(params.low_input)==4:
            self.B.minx.SetValue(str(params.low_input[0]))
            self.B.miny.SetValue(str(params.low_input[1]))
            self.B.minz.SetValue(str(params.low_input[2]))
            self.B.minr.SetValue(str(params.low_input[3]))
        if len(params.high_input)==4:
            self.B.maxx.SetValue(str(params.high_input[0]))
            self.B.maxy.SetValue(str(params.high_input[1]))
            self.B.maxz.SetValue(str(params.high_input[2]))
            self.B.maxr.SetValue(str(params.high_input[3]))

        #receptor data
        if params.receptor!="NA":
            self.R.usereceptor.SetValue(True)
            self.R.setstructure.SetValue(str(params.receptor))
        if len(params.low_input_rec)==4:
            self.R.minx.SetValue(str(params.low_input_rec[0]))
            self.R.miny.SetValue(str(params.low_input_rec[1]))
            self.R.minz.SetValue(str(params.low_input_rec[2]))
            self.R.minr.SetValue(str(params.low_input_rec[3]))
        if len(params.high_input_rec)==4:
            self.R.maxx.SetValue(str(params.high_input_rec[0]))
            self.R.maxy.SetValue(str(params.high_input_rec[1]))
            self.R.maxz.SetValue(str(params.high_input_rec[2]))
            self.R.maxr.SetValue(str(params.high_input_rec[3]))
        #self.add('z_padding','pad','int',10)


    def OnLaunch(self,event):
        """ Generate input file and call PSO"""
        r=self.check()
        if r == 0:
            self.makeFile(estimate=True)

            dia=LaunchDialog(self)
            dia.ShowModal()
            dia.Destroy()


    def OnVerify(self,event):
        """ Check that data in form is consistent with a PSO execution"""
        test=self.check()
        if test==0:
            self.errorPopup("Setup OK!")

    def check(self):

        if self.I.check()==-1:
            return -1

        if self.B.check()==-1:
            return -1

        if self.R.check()==-1:
            return -1

        if self.C.check()==-1:
            return -1

        if self.S.check()==-1:
            return -1

        return 0


    def errorPopup(self,msg):
        """ Popup error message"""
        dlg = wx.MessageDialog(self, msg, "Message!", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.


    def makeFile(self,fname=-1, estimate=False):
        """ read parameters and create input file """
        if fname==-1:
            fname='pso.txt'
        f=open(fname,'w')

        #write input structure data
        if str(self.I.editassemblystyle.GetValue())!="":
            f.write("style %s\n"%self.I.editassemblystyle.GetValue())
        if str(self.I.setdegree.GetValue())!="":
            f.write("degree %s\n"%self.I.setdegree.GetValue())

        if self.I.editassemblystyle.GetValue()=="flexible":
            if str(self.I.setstructure.GetValue())!="":
                f.write("trajectory %s\n"%self.I.setstructure.GetValue())
            if str(self.I.settop.GetValue())!="":
                f.write("topology %s\n"%self.I.settop.GetValue())
            if self.I.setproj.GetValue()!="":
                f.write("projection %s\n"%self.I.setproj.GetValue())
            if self.I.setratio.GetValue()!="":
                f.write("ratio %s\n"%self.I.setratio.GetValue())
        else:
            if self.I.setstructure.GetValue()!="":
                f.write("monomer %s\n"%self.I.setstructure.GetValue())

        #write constraint data
        if self.C.setconstraint.GetValue()!="":
            f.write("constraint %s\n"%self.C.setconstraint.GetValue())
        target=''
        keyword=''
        for x in xrange(0,self.C.lc.GetItemCount(),1):
            #min and max
            c=self.C.lc.GetItem(x, 1).GetText()
            if c=='':
                c="NaN"
            target+=c+" "
            #keyword
            k=self.C.lc.GetItem(x, 2).GetText()
            if k =='':
                k='NaN'
            else:
                keyword+=k+" "

        if self.C.lc.GetItemCount()>0:
            f.write("target %s\n"%(target))
            f.write("#keyword %s\n"%(keyword))

        if self.C.editmix.GetValue()!="":
            f.write("mixingWeight %s\n"%self.C.editmix.GetValue())

        #write boundary conditions
        if self.B.minx.GetValue()!="":
            minx=self.B.minx.GetValue()
        else:
            minx="NaN"
        if self.B.miny.GetValue()!="":
            miny=self.B.miny.GetValue()
        else:
            miny="NaN"
        if self.B.minz.GetValue()!="":
            minz=self.B.minz.GetValue()
        else:
            minz="NaN"
        if self.B.minr.GetValue()!="":
            minr=self.B.minr.GetValue()
        else:
            minr="NaN"
        if self.B.maxx.GetValue()!="":
            maxx=self.B.maxx.GetValue()
        else:
            maxx="NaN"
        if self.B.maxy.GetValue()!="":
            maxy=self.B.maxy.GetValue()
        else:
            maxy="NaN"
        if self.B.maxz.GetValue()!="":
            maxz=self.B.maxz.GetValue()
        else:
            maxz="NaN"
        if self.B.maxr.GetValue()!="":
            maxr=self.B.minr.GetValue()
        else:
            maxr="NaN"

        #write boundary data
        f.write("boundaryMin %s %s %s %s\n"%(minx,miny,minz,minr))
        f.write("boundaryMax %s %s %s %s\n"%(maxx,maxy,maxz,maxr))

        #write receptor conditions, if needed
        if self.R.usereceptor.Value==True:
            if self.R.minx.GetValue()!="":
                minx=self.R.minx.GetValue()
            else:
                minx="NaN"
            if self.R.miny.GetValue()!="":
                miny=self.R.miny.GetValue()
            else:
                miny="NaN"
            if self.B.minz.GetValue()!="":
                minz=self.R.minz.GetValue()
            else:
                minz="NaN"
            if self.R.minr.GetValue()!="":
                minr=self.R.minr.GetValue()
            else:
                minr="NaN"
            if self.R.maxx.GetValue()!="":
                maxx=self.R.maxx.GetValue()
            else:
                maxx="NaN"
            if self.R.maxy.GetValue()!="":
                maxy=self.R.maxy.GetValue()
            else:
                maxy="NaN"
            if self.R.maxz.GetValue()!="":
                maxz=self.R.maxz.GetValue()
            else:
                maxz="NaN"
            if self.R.maxr.GetValue()!="":
                maxr=self.R.minr.GetValue()
            else:
                maxr="NaN"

            f.write("boundaryMinReceptor %s %s %s %s\n"%(minx,miny,minz,minr))
            f.write("boundaryMaxReceptor %s %s %s %s\n"%(maxx,maxy,maxz,maxr))


        #write PSO setup data
        if str(self.S.editinertiamin)!='':
            f.write("inertiaMin %s\n"%self.S.editinertiamin)
        if str(self.S.editinertiamax)!='':
            f.write("inertiaMax %s\n"%self.S.editinertiamax)
        if str(self.S.editglobal)!='':
            f.write("cn %s\n"%self.S.editglobal)
        if str(self.S.editpersonal)!='':
            f.write("cp %s\n"%self.S.editpersonal)
        if str(self.S.editrep.GetValue())!='':
            f.write("repeat %s\n"%self.S.editrep.GetValue())
        if str(self.S.editsteps.GetValue())!='':
            f.write("steps %s\n"%self.S.editsteps.GetValue())
        if str(self.S.editpart.GetValue())!='':
            f.write("particles %s\n"%self.S.editpart.GetValue())
        if str(self.S.editneighsize)!='':
            f.write("neighborSize %s\n"%self.S.editneighsize)
        if str(self.S.editneigh)!='':
            f.write("neighborType %s\n"%self.S.editneigh)

        f.close()


class LaunchDialog(wx.Dialog):
    """ Ask the user hom many CPU to use """
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent,title="Launching module",size=(300,100))

        self.parent=parent

        self.nb=self.detectCPU()
        #selector
        self.text = wx.StaticText(self, label="Pick how many CPU to use: ")

        self.cpu=[]
        for x in xrange (1,min(int(self.parent.S.editpart.GetValue()), self.nb)+1,1):
            self.cpu.append(str(x))
        self.editcpu = wx.ComboBox(self, value=self.cpu[0], size=(50, -1), choices=self.cpu, style=wx.CB_DROPDOWN)

        boxpick = wx.BoxSizer(wx.HORIZONTAL)
        boxpick.Add(self.text, 0, wx.ALIGN_CENTER, 3 )
        boxpick.Add(self.editcpu, 1, wx.ALIGN_CENTER, 3 )

        #buttons
        self.ok =wx.Button(self, label="OK!")
        self.cancel =wx.Button(self, label="Cancel")

        boxbuttons = wx.BoxSizer(wx.HORIZONTAL)
        boxbuttons.Add(self.ok, 0, wx.ALIGN_CENTER, 3 )
        boxbuttons.Add(self.cancel, 1, wx.ALIGN_CENTER, 3 )

        #whole arrangement in the dialog
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(boxpick, 1, wx.ALIGN_CENTER, 3 )
        box.Add(boxbuttons, 1, wx.ALIGN_CENTER, 3 )
        self.SetSizer(box)
        self.Centre()

        #bindings
        self.Bind(wx.EVT_BUTTON, self.OnOk,self.ok,id=wx.ID_OK)
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)


    def detectCPU(self):
        """ detect the number of available CPU """
        if hasattr(os, "sysconf"):
            if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
                # Linux & Unix
                ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
                if isinstance(ncpus, int) and ncpus > 0:
                    return ncpus
            else:
                # OSX
                return int(os.popen2("sysctl -n hw.ncpu")[1].read())

        # Windows
        if os.environ.has_key("NUMBER_OF_PROCESSORS"):
            ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
            if ncpus > 0:
                return ncpus
        return 1

    def OnOk(self,e):
        pso_path=os.path.realpath(os.path.dirname(sys.argv[0]))
        print "> launching optimization with "+str(self.editcpu.GetValue())+" processors..."
        call='mpirun -n '+str(self.editcpu.GetValue())+' '+pso_path+'/POW.py DockSymmCircle pso.txt'
        subprocess.check_call(call,shell=True)
        os.remove('pso.txt')

        self.Close(True)

    def OnCancel(self,e):
        self.Close(True)


class IOPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        self.dirname=''
        self.parent=parent

        #title
        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="STRUCTURE")
        self.title.SetFont(font1)

        # neighbor type
        self.assemblystyle = ['rigid', 'flexible']
        self.lblassemblystyle = wx.StaticText(self, label="style:")
        self.editassemblystyle = wx.ComboBox(self, value=self.assemblystyle[0], size=(70, -1), choices=self.assemblystyle, style=wx.CB_DROPDOWN)

        #stoichiometry
        self.lbldegree=wx.StaticText(self, label="stoichiometry:")
        self.setdegree=wx.TextCtrl(self, value="",size=(35, -1))

        gs1 = wx.GridSizer(1, 4, 3, 3)
        gs1.AddMany([(self.lblassemblystyle, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editassemblystyle, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lbldegree, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.setdegree, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        #structure line
        self.lblstructure=wx.StaticText(self, label="structure:")
        self.setstructure=wx.TextCtrl(self, value="")
        self.structurebutton=wx.Button(self, label="select")

        #topology line
        self.lbltop=wx.StaticText(self, label="topology:")
        self.settop=wx.TextCtrl(self, value="")
        self.topbutton=wx.Button(self, label="select")
        self.settop.Enable(False)
        self.topbutton.Enable(False)

        #projection file
        self.lblproj=wx.StaticText(self, label="projection:")
        self.setproj=wx.TextCtrl(self, value="")
        self.projbutton=wx.Button(self, label="select")
        self.setproj.Enable(False)
        self.projbutton.Enable(False)

        #ratio
        self.lblratio=wx.StaticText(self, label="ratio:")
        self.setratio=wx.TextCtrl(self, value="0.9",size=(35, -1))
        self.setratio.Enable(False)

        #put structure, topology and projection stuff in a grid
        gs2 = wx.FlexGridSizer(3, 3, 3, 10)
        gs2.SetFlexibleDirection(wx.HORIZONTAL)
        gs2.AddGrowableCol(1,1)
        gs2.AddMany([(self.lblstructure, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setstructure, 1,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.structurebutton, 0, wx.ALIGN_LEFT, 3),
                    (self.lbltop, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3),
                    (self.settop, 1, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.topbutton, 0, wx.ALIGN_LEFT, 3),
                    (self.lblproj, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3),
                    (self.setproj, 1, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.projbutton, 0, wx.ALIGN_LEFT, 3),
                    (self.lblratio, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setratio, 0, wx.ALIGN_LEFT, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3)])

        self.Bind(wx.EVT_BUTTON, self.OnClickSelectStruct,self.structurebutton)
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectTop,self.topbutton)
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectProj,self.projbutton)
        self.Bind(wx.EVT_COMBOBOX, self.OnStyleCombobox,self.editassemblystyle)

        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs1, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs2, 1, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()

    def OnStyleCombobox(self,event):
        #if the style is rigid, disable everything, enable otherwise
        if self.editassemblystyle.Value=="rigid":
            self.activator=False
            #self.settop.Hide()
            #self.topbutton.Hide()
            #self.lbltop.Hide()
            #self.setproj.Hide()
            #self.projbutton.Hide()
            #self.lblproj.Hide()
            #self.setratio.Hide()
            #self.lblratio.Hide()
        else:
            self.activator=True
            #self.settop.Show()
            #self.topbutton.Show()
            #self.lbltop.Show()
            #self.setproj.Show()
            #self.projbutton.Show()
            #self.lblproj.Show()
            #self.setratio.Show()
            #self.lblratio.Show()

        self.settop.Enable(self.activator)
        self.topbutton.Enable(self.activator)
        self.setproj.Enable(self.activator)
        self.projbutton.Enable(self.activator)
        self.setratio.Enable(self.activator)

    def OnClickSelectStruct(self,event):
        dlg = wx.FileDialog(self, "Choose your structure (pdb or trajectory)", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setstructure.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()

    def OnClickSelectTop(self,event):
        dlg = wx.FileDialog(self, "Choose your topology file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.settop.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()

    def OnClickSelectProj(self,event):
        dlg = wx.FileDialog(self, "Choose your projection file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setproj.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


    def check(self):
        """ Verify parameters consistency"""

        #check stoichiometry value
        try:
            setdegree_f=int(self.setdegree.GetValue())
        except ValueError:
            self.parent.errorPopup("soichiometry should be an integer number greater than 1")
            return -1

        if setdegree_f<2:
            self.parent.errorPopup("assembly soichiometry should be greater or equal to 2")
            return -1

        #check constraint file existence
        fitfile=self.setstructure.GetValue()
        if fitfile == '':
            self.parent.errorPopup("structure not defined!")
            return -1
        if os.path.isfile(fitfile)==False:
            self.parent.errorPopup("structure does not exist!")
            return -1

        if self.editassemblystyle.GetValue()=="rigid":
            return 0

        #check topology file existence
        fitfile=self.settop.GetValue()
        if fitfile == '':
            self.parent.errorPopup("when performing flexible assembly, a topology file should be provided (crd or dcd format)!")
            return -1
        if os.path.isfile(fitfile)==False:
            self.parent.errorPopup("topology file does not exist!")
            return -1

        #check projection file existence
        projfile=self.setproj.GetValue()
        if projfile!='' and os.path.isfile(projfile)==False:
            self.parent.errorPopup("projections file does not exist!")
            return -1

        if projfile=='':
            try:
                setratio_f=float(self.setratio.GetValue())
            except ValueError:
                self.parent.errorPopup("ratio value undefined! Expecting a number within 0 and 1.")
                return -1
            if setratio_f<0 or setratio_f>1:
                self.parent.errorPopup("ratio should be a number within 0 and 1!")
                return -1

class ReceptorPanel(wx.Panel):
    def __init__(self, parent):

        self.dirname=''
        self.parent=parent

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.t = wx.StaticText(self, label="RECEPTOR")
        self.t.SetFont(font1)
        # use repulsion field
        self.usereceptor = wx.CheckBox(self, label="")

        self.title = wx.FlexGridSizer(1, 2, 3, 5)
        self.title.SetFlexibleDirection(wx.HORIZONTAL)
        self.title.AddMany([(self.t, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT,3),
                    (self.usereceptor, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 3)])

        #structure line
        self.lblstructure=wx.StaticText(self, label="structure:")
        self.setstructure=wx.TextCtrl(self, value="")
        self.structurebutton=wx.Button(self, label="select")

        #box the group
        gs1 = wx.FlexGridSizer(1, 3, 3, 5)
        gs1.SetFlexibleDirection(wx.HORIZONTAL)
        gs1.AddGrowableCol(1,1)
        gs1.AddMany([(self.lblstructure, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setstructure, 1,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.structurebutton, 0, wx.ALIGN_LEFT, 3)])

        self.Bind(wx.EVT_BUTTON, self.OnClickSelectReceptor,self.structurebutton)

        #top lables
        self.lblx = wx.StaticText(self, label="rot. x")
        self.lbly = wx.StaticText(self, label="rot. y")
        self.lblz = wx.StaticText(self, label="rot. z")
        self.lblr = wx.StaticText(self, label="trans. z")

        #min boundaries
        self.lblmin = wx.StaticText(self, label="min :")
        self.minx = wx.TextCtrl(self, value="0", size=(40,-1))
        self.miny = wx.TextCtrl(self, value="0", size=(40,-1))
        self.minz = wx.TextCtrl(self, value="0", size=(40,-1))
        self.minr = wx.TextCtrl(self, value="", size=(40,-1))

        #max boundaries
        self.lblmax = wx.StaticText(self, label="max :")
        self.maxx = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxy = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxz = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxr = wx.TextCtrl(self, value="", size=(40,-1))

        #grid of PSO constants
        gs = wx.GridSizer(3, 5, 3, 3)
        gs.AddMany([(wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lblx, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lbly, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblz, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblr, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblmin, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minx, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.miny, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minz, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minr, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblmax, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxx, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxy, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxz, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxr, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        self.Bind(wx.EVT_CHECKBOX, self.OnUseReceptor,self.usereceptor)
        self.structurebutton.Enable(False)
        self.setstructure.Enable(False)
        self.minx.Enable(False)
        self.miny.Enable(False)
        self.minz.Enable(False)
        self.minr.Enable(False)
        self.maxx.Enable(False)
        self.maxy.Enable(False)
        self.maxz.Enable(False)
        self.maxr.Enable(False)


        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        #box.Add(self.usereceptor, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs1, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        #box.Add(self.repel, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()

    def OnClickSelectReceptor(self,event):
        dlg = wx.FileDialog(self, "Choose your receptor file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setstructure.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


    def OnUseReceptor(self,event):
        #if the style is rigid, disable everything, enable otherwise
        if self.usereceptor.Value==False:
            self.activator=False
        else:
            self.activator=True

        self.structurebutton.Enable(self.activator)
        self.setstructure.Enable(self.activator)
        self.minx.Enable(self.activator)
        self.miny.Enable(self.activator)
        self.minz.Enable(self.activator)
        self.minr.Enable(self.activator)
        self.maxx.Enable(self.activator)
        self.maxy.Enable(self.activator)
        self.maxz.Enable(self.activator)
        self.maxr.Enable(self.activator)


    def check(self):
        """ Verify parameters consistency"""

        if self.usereceptor.GetValue()==False:
            return 0

        #check constraint file existence
        self.setstructure.GetValue()
        fitfile=self.setstructure.GetValue()
        if fitfile == '':
            self.parent.errorPopup("receptor structure not defined!")
            return -1
        if os.path.isfile(fitfile)==False:
            self.parent.errorPopup("receptor structure does not exist!")
            return -1

        #type consistency for min and max boundaries
        try:
            minx_f=float(self.minx.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor x min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            miny_f=float(self.miny.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor y min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            minz_f=float(self.minz.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor z min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            minr_f=float(self.minr.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor radius min boundary must be defined as a real number")
            return -1

        try:
            maxx_f=float(self.maxx.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor x min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            maxy_f=float(self.maxy.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor y max boundary must be defined as a number within 0 and 360")
            return -1

        try:
            maxz_f=float(self.maxz.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor z max boundary must be defined as a number within 0 and 360")
            return -1

        try:
            maxr_f=float(self.maxr.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor radius max boundary must be defined as a real number")
            return -1

        #value consistency for min and max boundaries
        if minx_f>maxx_f:
            self.parent.errorPopup("receptor in rotation around x axis, min boundary is greater than max!")
            return -1

        if miny_f>maxy_f:
            self.parent.errorPopup("receptor in rotation around y axis, min boundary is greater than max!")
            return -1

        if minz_f>maxz_f:
            self.parent.errorPopup("receptor in rotation around z axis, min boundary is greater than max!")
            return -1

        if minr_f>maxr_f:
            self.parent.errorPopup("receptor in assembly radius definition, min boundary is greater than max!")
            return -1


class ConstraintPanel(wx.Panel):
    def __init__(self, parent):

        self.dirname=''
        self.parent=parent

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        self.font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)

        self.title = wx.StaticText(self, label="FITNESS")
        self.title.SetFont(self.font1)

        #fitness line
        self.lblconstraint=wx.StaticText(self, label="constraint:")
        self.setconstraint=wx.TextCtrl(self, value="")
        self.constraintbutton=wx.Button(self, label="select")
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectConstraint,self.constraintbutton)

        # boundary conditions
        self.lc = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
        self.lc.InsertColumn(0, 'ID')
        self.lc.InsertColumn(1, 'value')
        self.lc.InsertColumn(2, 'keyword')
        self.lc.SetColumnWidth(0, 50)

        #boundary conditions buttons box
        self.buttonsbox = wx.BoxSizer(wx.VERTICAL)
        self.buttonsbox.Add(wx.Button(self, 1, "Add"), 0 )
        self.buttonsbox.Add(wx.Button(self, 2, "Edit"), 0 )
        self.buttonsbox.Add(wx.Button(self, 3, "Remove"), 0 )

        #box the group
        gs = wx.FlexGridSizer(2, 3, 3, 5)
        #gs.SetFlexibleDirection(wx.HORIZONTAL)
        gs.AddGrowableCol(1,1)
        gs.AddGrowableRow(1,1)
        gs.AddMany([(self.lblconstraint, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setconstraint, 0, wx.ALIGN_TOP | wx.ALIGN_LEFT| wx.EXPAND, 3),
                    (self.constraintbutton, 0, wx.ALIGN_LEFT, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lc, 1, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.buttonsbox, 0, wx.ALIGN_TOP | wx.ALIGN_LEFT, 3)])

        # energy vs geometry mixing weight
        self.lblmix = wx.StaticText(self, label="mix coeff.:")
        self.editmix = wx.TextCtrl(self, value="0.2", size=(40,-1))

        #log filtering criteria
        self.lblfilter = wx.StaticText(self, label="filtering thresh.:")
        self.editfilter = wx.TextCtrl(self, value="0", size=(40,-1))

        gs2 = wx.GridSizer(1, 4, 3, 3)
        gs2.AddMany([(self.lblmix, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editmix, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblfilter, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editfilter, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        #draw title, boundaries and buttons list
        boundbox = wx.BoxSizer(wx.VERTICAL)
        boundbox.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        boundbox.Add(gs, 1, wx.EXPAND | wx.ALL, 3 )
        boundbox.Add(gs2, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(boundbox,1)
        self.SetSize((400,600))
        self.Centre()

        #buttons bindings
        self.Bind(wx.EVT_BUTTON, self.OnBadd,id=1)
        self.Bind(wx.EVT_BUTTON, self.OnBedit,id=2)
        self.Bind(wx.EVT_BUTTON, self.OnBremove,id=3)


    def OnClickSelectConstraint(self,event):
        dlg = wx.FileDialog(self, "Choose your constraint file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setconstraint.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


    def OnBadd(self,event):

        dia=Boundaries(self)
        dia.ShowModal()
        dia.Destroy()


    def OnBedit(self,event):

        index = self.lc.GetFocusedItem()

        params=[]
        params.append(index)
        if index != -1 :
            params.append(self.lc.GetItem(index, 1).GetText())
            params.append(self.lc.GetItem(index, 2).GetText())
            params.append(self.lc.GetItem(index, 3).GetText())
            params.append(self.lc.GetItem(index, 4).GetText())
            params.append(self.lc.GetItem(index, 5).GetText())

            dia=Boundaries(self,params)
            dia.ShowModal()
            dia.Destroy()

    def OnBremove(self,event):
        index = self.lc.GetFocusedItem()
        self.lc.DeleteItem(index)

        #renumber items after deletion
        num_items = self.lc.GetItemCount()
        for x in xrange(index,num_items,1):
            self.lc.SetStringItem(x, 0, str(x+1))

    def check(self):
        #check constraint file existence
        self.setconstraint.GetValue()
        fitfile=self.setconstraint.GetValue()
        if fitfile == '':
            self.parent.errorPopup("constraint file not defined!")
            return -1
        if os.path.isfile(fitfile)==False:
            self.parent.errorPopup("constraint file does not exist!")
            return -1

        l=self.lc.GetItemCount()
        if l == 0:
            self.parent.errorPopup("no constraints have been defined!")
            return -1
        #test that the number of values outputted by the function and the number of lines is the same

        for x in xrange(0,l,1):
            #check min and max
            m1=self.lc.GetItem(x, 1).GetText()
            if m1=='':
                self.parent.errorPopup(" constraint %s has no value set!"%(x+1))
                return -1
            try:
                m1_f=float(m1)
            except ValueError:
                self.parent.errorPopup("constraint %s should be a number"%(x+1))
                return -1

        try:
            mix_f=float(self.editmix.GetValue())
        except ValueError:
            self.parent.errorPopup("mixing coefficient should be a number!")
            return -1
        if mix_f<0 or mix_f>1:
            self.parent.errorPopup("mixing coefficient should be a number within 0 and 1")
            return -1


class BoundaryPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        self.parent=parent

        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="BOUNDARIES")
        self.title.SetFont(font1)

        #top lables
        self.lblx = wx.StaticText(self, label="rot. x")
        self.lbly = wx.StaticText(self, label="rot. y")
        self.lblz = wx.StaticText(self, label="rot. z")
        self.lblr = wx.StaticText(self, label="radius")

        #min boundaries
        self.lblmin = wx.StaticText(self, label="min :")
        self.minx = wx.TextCtrl(self, value="0", size=(40,-1))
        self.miny = wx.TextCtrl(self, value="0", size=(40,-1))
        self.minz = wx.TextCtrl(self, value="0", size=(40,-1))
        self.minr = wx.TextCtrl(self, value="", size=(40,-1))

        #max boundaries
        self.lblmax = wx.StaticText(self, label="max :")
        self.maxx = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxy = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxz = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxr = wx.TextCtrl(self, value="", size=(40,-1))

        #grid of PSO constants
        gs = wx.GridSizer(3, 5, 3, 3)
        gs.AddMany([(wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lblx, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lbly, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblz, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblr, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblmin, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minx, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.miny, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minz, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minr, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblmax, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxx, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxy, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxz, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxr, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        #box.Add(self.repel, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box,1)
        self.Centre()

    def check(self):
        """ Verify parameters consistency"""

        #type consistency for min and max boundaries
        try:
            minx_f=float(self.minx.GetValue())
        except ValueError:
            self.parent.errorPopup("x min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            miny_f=float(self.miny.GetValue())
        except ValueError:
            self.parent.errorPopup("y min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            minz_f=float(self.minz.GetValue() )
        except ValueError:
            self.parent.errorPopup("z min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            minr_f=float(self.minr.GetValue())
        except ValueError:
            self.parent.errorPopup("radius min boundary must be defined as a real number")
            return -1

        try:
            maxx_f=float(self.maxx.GetValue())
        except ValueError:
            self.parent.errorPopup("x min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            maxy_f=float(self.maxy.GetValue())
        except ValueError:
            self.parent.errorPopup("y max boundary must be defined as a number within 0 and 360")
            return -1

        try:
            maxz_f=float(self.maxz.GetValue())
        except ValueError:
            self.parent.errorPopup("z max boundary must be defined as a number within 0 and 360")
            return -1

        try:
            maxr_f=float(self.maxr.GetValue())
        except ValueError:
            self.parent.errorPopup("radius max boundary must be defined as a real number")
            return -1

        #value consistency for min and max boundaries
        if minx_f>maxx_f:
            self.parent.errorPopup("in rotation around x axis, min boundary is greater than max!")
            return -1

        if miny_f>maxy_f:
            self.parent.errorPopup("in rotation around y axis, min boundary is greater than max!")
            return -1

        if minz_f>maxz_f:
            self.parent.errorPopup("in rotation around z axis, min boundary is greater than max!")
            return -1

        if minr_f>maxr_f:
            self.parent.errorPopup("in assembly radius definition, min boundary is greater than max!")
            return -1



class ProtocolPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        self.parent=parent

        #PSO behavior values
        self.editneighsize=1
        self.editneigh="geographic"
        self.editinertiamin=0.3
        self.editinertiamax=0.7
        self.editpersonal=1.2
        self.editglobal=1.4

        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="PROTOCOL")
        self.title.SetFont(font1)

        # PSO steps
        self.lblsteps = wx.StaticText(self, label="Steps")
        self.editsteps = wx.TextCtrl(self, value="100", size=(40,-1))
        # particles
        self.lblpart = wx.StaticText(self, label="Particles")
        self.editpart = wx.TextCtrl(self, value="40", size=(40,-1))
        #protocol
        self.lblrep = wx.StaticText(self, label="Repetitions")
        self.editrep = wx.TextCtrl(self, value="1", size=(40,-1))

        self.detailsbutton=wx.Button(self, label="Details...")
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectDetails,self.detailsbutton)

        gs = wx.GridSizer(2, 4, 3, 3)
        gs.AddMany([(self.lblsteps, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblpart, 0,  wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblrep, 0,  wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.editsteps, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editpart, 0,  wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editrep, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.detailsbutton, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3)])

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()

    def OnClickSelectDetails(self,event):
        params=[self.editneigh,self.editneighsize,self.editinertiamin,self.editinertiamax,self.editpersonal,self.editglobal]
        dia=PSO(self,params)
        dia.ShowModal()
        dia.Destroy()


    def check(self):
        #check inertia
        imin=self.editinertiamin
        imax=self.editinertiamax
        try:
            imin_f=float(imin)
        except ValueError:
            self.parent.errorPopup("min inertia must be defined as a number within 0 and 1")
            return -1
        try:
            imax_f=float(imax)
        except ValueError:
            self.parent.errorPopup("max inertia must be defined as a number within 0 and 1")
            return -1
        if imin_f<0:
            self.parent.errorPopup("min inertia should be greater than zero")
            return -1
        if imax_f>=1:
            self.parent.errorPopup("max inertia should be smaller than 1")
            return -1
        if imin_f>imax_f:
            self.parent.errorPopup("max inertia should be greater or equal than min inertia")
            return -1

        cn=self.editglobal
        try:
            cn_f=float(cn)
        except ValueError:
            self.parent.errorPopup("global weight should be a number")
            return -1
        if cn_f<0:
            self.parent.errorPopup("global weight should be positive")
            return -1

        cp=self.editpersonal
        try:
            cp_f=float(cp)
        except ValueError:
            self.parent.errorPopup("personal weight should be a number")
            return -1
        if cp_f<0:
            self.parent.errorPopup("personal weight should be positive")
            return -1

        psosteps=self.editsteps.GetValue()
        try:
            psosteps_f=int(psosteps)
        except ValueError:
            self.parent.errorPopup("number of PSO steps should be an integer")
            return -1
        if psosteps_f<0:
            self.parent.errorPopup("number of PSO steps should be greater than zero")
            return -1

        rep=self.editrep.GetValue()
        try:
            rep_f=int(rep)
        except ValueError:
            self.parent.errorPopup("number of PSO repetitions should be an integer")
            return -1
        if rep_f<0:
            self.parent.errorPopup("number of PSO repetitions should be greater than zero")
            return -1

        particles=self.editpart.GetValue()
        try:
            particles_f=int(particles)
        except ValueError:
            self.parent.errorPopup("number of particles should be an integer")
            return -1
        if particles_f<0:
            self.parent.errorPopup("number of particles should be greater than zero")
            return -1

        nsize=self.editneighsize
        try:
            nsize_f=int(nsize)
        except ValueError:
            self.parent.errorPopup("neighborhood size should be an integer")
            return -1
        if nsize_f<0:
            self.parent.errorPopup("neighborhood size should be greater than zero")
            return -1
        if nsize_f>=particles_f:
            self.parent.errorPopup("neighborhood size should be smaller than number of particles")
            return -1

        ntype=self.editneigh
        if ntype!='geographic' and ntype!='indexed':
            self.parent.errorPopup("neighborhood type should be either geographic, either indexed")
            return -1



class PSO(wx.Dialog):
    def __init__(self, parent,p):

        wx.Dialog.__init__(self, parent,title="PSO parametrization",size=(350,125))

        self.parent=parent

        # neighbor type
        self.neighList = ['geographic', 'indexed']
        self.lblneigh = wx.StaticText(self, label="Neighbor style : ")
        self.editneigh = wx.ComboBox(self, value=str(p[0]), size=(80, -1), choices=self.neighList, style=wx.CB_DROPDOWN)
        # neighborhood size
        self.lblneighsize = wx.StaticText(self, label="Neighbors nb. :")
        self.editneighsize = wx.TextCtrl(self, value=str(p[1]), size=(40,-1))

        self.lblinertiamin = wx.StaticText(self, label="inertia min :")
        self.editinertiamin = wx.TextCtrl(self, value=str(p[2]), size=(40,-1))
        self.lblinertiamax = wx.StaticText(self, label="inertia max :")
        self.editinertiamax = wx.TextCtrl(self, value=str(p[3]), size=(40,-1))

        self.lblpersonal = wx.StaticText(self, label="personal weight :")
        self.editpersonal = wx.TextCtrl(self, value=str(p[4]), size=(40,-1))
        self.lblglobal = wx.StaticText(self, label="global weight :")
        self.editglobal = wx.TextCtrl(self, value=str(p[5]), size=(40,-1))

        # use repulsion field
        #self.repel = wx.CheckBox(self, label="use repulsion field?")

        #grid of PSO constants
        gs = wx.GridSizer(3, 4, 3, 3)
        gs.AddMany([(self.lblneigh, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editneigh, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblneighsize, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editneighsize, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblinertiamin, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editinertiamin, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblinertiamax, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editinertiamax, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblpersonal, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editpersonal, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblglobal, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editglobal, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])


        self.okbutton=wx.Button(self, label="OK!")
        self.Bind(wx.EVT_BUTTON, self.OnOk,self.okbutton,id=wx.ID_OK)
        self.cancel =wx.Button(self, label="Cancel")
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)
        self.buttonbox = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonbox.Add(self.okbutton, 0, wx.RIGHT, 3 )
        self.buttonbox.Add(self.cancel, 0, wx.LEFT, 3 )

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.buttonbox, 1, wx.CENTER, 3 )
        self.SetSizer(box)
        self.Centre()

    def OnOk(self,e):
        self.parent.editneighsize=self.editneighsize.GetValue()
        self.parent.editneigh=self.editneigh.GetValue()
        self.parent.editinertiamin=self.editinertiamin.GetValue()
        self.parent.editinertiamax=self.editinertiamax.GetValue()
        self.parent.editpersonal=self.editpersonal.GetValue()
        self.parent.editglobal=self.editglobal.GetValue()

        self.Close(True)  # Close the frame.

    def OnCancel(self,e):
        self.Close(True)


class Boundaries(wx.Dialog):
    def __init__(self, parent,params=-1):
        wx.Dialog.__init__(self, parent,title="Boundaries Editor",size=(300,75))

        self.params=params
        self.parent=parent

        self.lblval = wx.StaticText(self, label="value:")
        self.editval = wx.TextCtrl(self, value="", size=(40,-1))
        self.lblkeyword = wx.StaticText(self, label="keyword:")
        self.editkeyword = wx.TextCtrl(self, value="", size=(75,-1))

        if params!=-1:
            self.editval.SetValue(self.params[1])
            self.editkeyword.SetValue(self.params[2])

        gs = wx.GridSizer(1, 4, 3, 3)
        gs.AddMany([(self.lblval, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editval, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblkeyword, 1,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editkeyword, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        self.okbutton=wx.Button(self, label="OK!")
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.okbutton,id=wx.ID_OK)
        self.cancel =wx.Button(self, label="Cancel")
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)
        self.buttonbox = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonbox.Add(self.okbutton, 0, wx.RIGHT, 3 )
        self.buttonbox.Add(self.cancel, 0, wx.LEFT, 3 )

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs, 0, wx.EXPAND, 3 )
        box.Add(self.buttonbox, 1, wx.CENTER, 3 )
        self.SetSizer(box)
        self.Centre()

    def OnOk(self,e):

        if self.params!=-1:
            num_items=self.params[0]
        else:
            num_items = self.parent.lc.GetItemCount()
            self.parent.lc.InsertStringItem(num_items, str(num_items+1))

        self.parent.lc.SetStringItem(num_items, 1, self.editval.GetValue())
        self.parent.lc.SetStringItem(num_items, 2, self.editkeyword.GetValue())

        self.Close(True)  # Close the frame.

    def OnCancel(self,e):
        self.Close(True)  # Close the frame.


app = wx.App(False)
frame = MainWindow(None, "POW Assembly of Symmetrical Assemblies")
frame.Show()
app.MainLoop()
#! /usr/bin/env python
# Copyright (c) 2012 EPFL (Ecole Polytechnique federale de Lausanne)
# Laboratory for Biomolecular Modeling, School of Life Sciences
#
# POW is free software ;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation ;
# either version 2 of the License, or (at your option) any later version.
# POW is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with POW ;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
#
# Author : Matteo Degiacomi, matteothomas.degiacomi@epfl.ch
# Web site : http://lbm.epfl.ch


import os, sys
import wx
import subprocess
import DockSymmCircle as module

class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        self.dirname=''
        # A "-1" in the size parameter instructs wxWidgets to use the default size.
        # In this case, we select 200px width and the default height.
        wx.Frame.__init__(self, parent, title=title, size=(600,500))
        self.CreateStatusBar() # A Status bar in the bottom of the window


        self.SetBackgroundColour("white")
        
        # Setting up the menu bar
        filemenu= wx.Menu()
        menuLoad = filemenu.Append(wx.ID_OPEN, "&Load..."," Load a setup file")
        menuSave = filemenu.Append(wx.ID_SAVEAS, "&Save as..."," Save the current setup")
        menuHelp = filemenu.Append(wx.ID_HELP, "&Help"," Optimization user manual")
        menuAbout= filemenu.Append(wx.ID_ABOUT, "&About"," Information about this program")
        filemenu.AppendSeparator()
        menuExit = filemenu.Append(wx.ID_EXIT,"E&xit"," Terminate the program")

        # Creating the menu bar
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

        # launch and verify buttons    
        self.launch=wx.Button(self, label="launch!")
        self.verify=wx.Button(self, label="verify")
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        buttons.Add(self.verify, 0 ) 
        buttons.Add(self.launch, 0 )    

        #panels        
        self.I=IOPanel(self)      
        self.C=ConstraintPanel(self)
        self.B=BoundaryPanel(self)
        self.R=ReceptorPanel(self)  
        self.S=ProtocolPanel(self)

        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.I, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.B, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.R, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.C, 1, wx.EXPAND | wx.ALL, 3 )
        #box.Add(self.P, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.S, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(buttons, 0, wx.ALIGN_CENTER)
        self.SetSizer(box)
        self.SetSize((400,800))
        self.Centre()  

        # bind events to menu items
        self.Bind(wx.EVT_MENU, self.OnLoad, menuLoad)
        self.Bind(wx.EVT_MENU, self.OnSave, menuSave)
        self.Bind(wx.EVT_MENU, self.OnHelp, menuHelp)
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
        self.Bind(wx.EVT_BUTTON, self.OnLaunch,self.launch)
        self.Bind(wx.EVT_BUTTON, self.OnVerify,self.verify)
        
        self.Show()


    def OnAbout(self,e):
        # Create a message dialog box
        dlg = wx.MessageDialog(self, " Prediction of circular protein assemblies@POW \n Matteo Degiacomi \n 2012", "About POW", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.

    def OnExit(self,e):
        """ Close the interface"""
        self.Close(True)  # Close the frame.

    def OnSave(self,e):
        """ Save a file"""
        dlg = wx.FileDialog(self, "Chose a name for your setup file!", self.dirname, "", "*.*", wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.makeFile(os.path.join(self.dirname, self.filename))
        dlg.Destroy()

    def OnHelp(self,e):
        """ Open lbm webpage """
        dlg = wx.MessageDialog(self, "Code updates and user manual available at \n lbm.epf.ch/resources", "POW help", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.


    def OnLoad(self,e):
        """ Open a file"""
        infile=-1
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            infile =os.path.join(self.dirname, self.filename)
        dlg.Destroy()
        
        #call file parsing
        if infile!=-1:
            if os.path.isfile(infile):
                self.parse(infile)
        
    def init(self):
        """ Reset all the form fields"""
        self.S.editsteps.SetValue('')
        self.S.editpart.SetValue('')
        self.S.editrep.SetValue('')
        self.C.setconstraint.SetValue('')        
        for __ in xrange(0,self.C.lc.GetItemCount(),1):
            self.C.lc.DeleteItem()
        
    def parse(self,infile):
        """ Input file parser"""
        
        self.init()
        
        params=module.Parser() #read user defined variables
	params.add_standard() #add default variables
	params.set_default_values() #set default values to all defined variables
	params.parse(infile) #parse input file
	
	#protocol panel
        self.S.editrep.SetValue(str(params.repeat))
        self.S.editsteps.SetValue(str(params.max_steps))
        self.S.editpart.SetValue(str(params.n_particles))

	#protocol advanced settings
        self.S.editneighsize=str(params.neigh_size)
        self.S.editneigh=str(params.neigh_type)
        self.S.editinertiamax=float(params.inertia_max)
        self.S.editinertiamin=float(params.inertia_min)
        self.S.editpersonal=float(params.cp)
        self.S.editglobal=float(params.cn)

	#constraint panel
        #insert data from target values list
        for x in xrange(0,len(params.target),1):
            self.C.lc.InsertStringItem(x, str(x+1))
            if str(params.target[x]) != "NaN":
                self.C.lc.SetStringItem(x, 1, str(params.target[x]))
            #if len(params.keyword)!=0 and str(params.keyword(x)) != "NaN":
            #    self.B.lc.SetStringItem(x, 2, str(params.keyword(x)))        

        self.C.setconstraint.SetValue(str(params.constraint))
        self.C.editmix.SetValue(str(params.mix_weight))
        self.C.editfilter.SetValue(str(params.accept))
        #self.add('detectClash','detect_clash','str',"on")

	#input panel
	self.I.setdegree.SetValue(str(params.degree))
	self.I.editassemblystyle.SetValue(str(params.style))
	self.I.setstructure.SetValue(str(params.pdb_file_name))
	self.I.settop.SetValue(str(params.topology))
	self.I.setproj.SetValue(str(params.proj_file))
	self.I.setratio.SetValue(str(params.ratio))
       	#self.add('align','align','str',"no")
	#self.add('trajectory','trajectory','str',"NA")
	#self.add('trajSelection','trajselection','str',"NA")

	print str(params.pdb_file_name)

	#boundary panel
	if len(params.low_input)==4:
		self.B.minx.SetValue(str(params.low_input[0]))
		self.B.miny.SetValue(str(params.low_input[1]))
		self.B.minz.SetValue(str(params.low_input[2]))
		self.B.minr.SetValue(str(params.low_input[3]))
	if len(params.high_input)==4:
		self.B.maxx.SetValue(str(params.high_input[0]))
		self.B.maxy.SetValue(str(params.high_input[1]))
		self.B.maxz.SetValue(str(params.high_input[2]))
		self.B.maxr.SetValue(str(params.high_input[3]))

	#receptor data
        if params.receptor!="NA":
		self.R.usereceptor.SetValue(True)
		self.R.setstructure.SetValue(str(params.receptor))
	if len(params.low_input_rec)==4:
		self.R.minx.SetValue(str(params.low_input_rec[0]))
		self.R.miny.SetValue(str(params.low_input_rec[1]))
		self.R.minz.SetValue(str(params.low_input_rec[2]))
		self.R.minr.SetValue(str(params.low_input_rec[3]))
	if len(params.high_input_rec)==4:
		self.R.maxx.SetValue(str(params.high_input_rec[0]))
		self.R.maxy.SetValue(str(params.high_input_rec[1]))
		self.R.maxz.SetValue(str(params.high_input_rec[2]))
		self.R.maxr.SetValue(str(params.high_input_rec[3]))
	#self.add('z_padding','pad','int',10)


    def OnLaunch(self,event):
        """ Generate input file and call PSO"""
        r=self.check()
        if r == 0:
            self.makeFile(estimate=True)
                        
            dia=LaunchDialog(self)
            dia.ShowModal()
            dia.Destroy()
        
        
    def OnVerify(self,event):
        """ Check that data in form is consistent with a PSO execution"""
        test=self.check()
        if test==0:
            self.errorPopup("Setup OK!")

    def check(self):
        
        if self.I.check()==-1:
            return -1
           
        if self.B.check()==-1:
            return -1
        
        if self.R.check()==-1:
            return -1
        
        if self.C.check()==-1:
            return -1        
        
        if self.S.check()==-1:
            return -1
                   
        return 0


    def errorPopup(self,msg):
        """ Popup error message"""
        dlg = wx.MessageDialog(self, msg, "Message!", wx.OK)
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.


    def makeFile(self,fname=-1, estimate=False):
        """ read parameters and create input file """
        if fname==-1:
            fname='pso.txt'            
        f=open(fname,'w')
        
        #write input structure data
        if str(self.I.editassemblystyle.GetValue())!="":
            f.write("style %s\n"%self.I.editassemblystyle.GetValue())
        if str(self.I.setdegree.GetValue())!="":
            f.write("degree %s\n"%self.I.setdegree.GetValue())
        
        if self.I.editassemblystyle.GetValue()=="flexible":
            if str(self.I.setstructure.GetValue())!="":            
                f.write("trajectory %s\n"%self.I.setstructure.GetValue())
            if str(self.I.settop.GetValue())!="":     
                f.write("topology %s\n"%self.I.settop.GetValue())
            if self.I.setproj.GetValue()!="":
                f.write("projection %s\n"%self.I.setproj.GetValue())     
            if self.I.setratio.GetValue()!="":
                f.write("ratio %s\n"%self.I.setratio.GetValue())
        else:
            if self.I.setstructure.GetValue()!="":
                f.write("monomer %s\n"%self.I.setstructure.GetValue())
        
        #write constraint data
        if self.C.setconstraint.GetValue()!="":
            f.write("constraint %s\n"%self.C.setconstraint.GetValue())
        target=''
        keyword=''
        for x in xrange(0,self.C.lc.GetItemCount(),1):
            #min and max
            c=self.C.lc.GetItem(x, 1).GetText()
            if c=='':
                c="NaN"
            target+=c+" "
            #keyword
            k=self.C.lc.GetItem(x, 2).GetText()
            if k =='':
                k='NaN'
            else:
                keyword+=k+" "
    
        if self.C.lc.GetItemCount()>0:
            f.write("target %s\n"%(target))        
            f.write("#keyword %s\n"%(keyword))        
        
        if self.C.editmix.GetValue()!="":
            f.write("mixingWeight %s\n"%self.C.editmix.GetValue())
        
        #write boundary conditions
        if self.B.minx.GetValue()!="":
            minx=self.B.minx.GetValue()
        else:
            minx="NaN"
        if self.B.miny.GetValue()!="":
            miny=self.B.miny.GetValue()
        else:
            miny="NaN"
        if self.B.minz.GetValue()!="":
            minz=self.B.minz.GetValue()
        else:
            minz="NaN"
        if self.B.minr.GetValue()!="":
            minr=self.B.minr.GetValue()
        else:
            minr="NaN"
        if self.B.maxx.GetValue()!="":
            maxx=self.B.maxx.GetValue()
        else:
            maxx="NaN"
        if self.B.maxy.GetValue()!="":
            maxy=self.B.maxy.GetValue()
        else:
            maxy="NaN"
        if self.B.maxz.GetValue()!="":
            maxz=self.B.maxz.GetValue()
        else:
            maxz="NaN"
        if self.B.maxr.GetValue()!="":
            maxr=self.B.minr.GetValue()
        else:
            maxr="NaN"
                     
        #write boundary data
        f.write("boundaryMin %s %s %s %s\n"%(minx,miny,minz,minr))
        f.write("boundaryMax %s %s %s %s\n"%(maxx,maxy,maxz,maxr))

        #write receptor conditions, if needed
        if self.R.usereceptor.Value==True:
            if self.R.minx.GetValue()!="":
                minx=self.R.minx.GetValue()
            else:
                minx="NaN"
            if self.R.miny.GetValue()!="":
                miny=self.R.miny.GetValue()
            else:
                miny="NaN"
            if self.B.minz.GetValue()!="":
                minz=self.R.minz.GetValue()
            else:
                minz="NaN"
            if self.R.minr.GetValue()!="":
                minr=self.R.minr.GetValue()
            else:
                minr="NaN"
            if self.R.maxx.GetValue()!="":
                maxx=self.R.maxx.GetValue()
            else:
                maxx="NaN"
            if self.R.maxy.GetValue()!="":
                maxy=self.R.maxy.GetValue()
            else:
                maxy="NaN"
            if self.R.maxz.GetValue()!="":
                maxz=self.R.maxz.GetValue()
            else:
                maxz="NaN"
            if self.R.maxr.GetValue()!="":
                maxr=self.R.minr.GetValue()
            else:
                maxr="NaN"
    
            f.write("boundaryMinReceptor %s %s %s %s\n"%(minx,miny,minz,minr))
            f.write("boundaryMaxReceptor %s %s %s %s\n"%(maxx,maxy,maxz,maxr))


        #write PSO setup data
        if str(self.S.editinertiamin)!='':
            f.write("inertiaMin %s\n"%self.S.editinertiamin)
        if str(self.S.editinertiamax)!='':
            f.write("inertiaMax %s\n"%self.S.editinertiamax)
        if str(self.S.editglobal)!='':
            f.write("cn %s\n"%self.S.editglobal)
        if str(self.S.editpersonal)!='':
            f.write("cp %s\n"%self.S.editpersonal)
        if str(self.S.editrep.GetValue())!='':
            f.write("repeat %s\n"%self.S.editrep.GetValue())
        if str(self.S.editsteps.GetValue())!='':
            f.write("steps %s\n"%self.S.editsteps.GetValue())
        if str(self.S.editpart.GetValue())!='':
            f.write("particles %s\n"%self.S.editpart.GetValue())
        if str(self.S.editneighsize)!='':
            f.write("neighborSize %s\n"%self.S.editneighsize)
        if str(self.S.editneigh)!='':
            f.write("neighborType %s\n"%self.S.editneigh)

        f.close()
        
        
class LaunchDialog(wx.Dialog):
    """ Ask the user hom many CPU to use """
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent,title="Launching module",size=(300,100))
  
        self.parent=parent
        
        self.nb=self.detectCPU()
        #selector
        self.text = wx.StaticText(self, label="Pick how many CPU to use: ")

        self.cpu=[]
        for x in xrange (1,min(int(self.parent.S.editpart.GetValue()), self.nb)+1,1):            
            self.cpu.append(str(x))
        self.editcpu = wx.ComboBox(self, value=self.cpu[0], size=(50, -1), choices=self.cpu, style=wx.CB_DROPDOWN)
        
        boxpick = wx.BoxSizer(wx.HORIZONTAL)
        boxpick.Add(self.text, 0, wx.ALIGN_CENTER, 3 )
        boxpick.Add(self.editcpu, 1, wx.ALIGN_CENTER, 3 )

        #buttons
        self.ok =wx.Button(self, label="OK!")
        self.cancel =wx.Button(self, label="Cancel")
        
        boxbuttons = wx.BoxSizer(wx.HORIZONTAL)
        boxbuttons.Add(self.ok, 0, wx.ALIGN_CENTER, 3 )
        boxbuttons.Add(self.cancel, 1, wx.ALIGN_CENTER, 3 )

        #whole arrangement in the dialog
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(boxpick, 1, wx.ALIGN_CENTER, 3 )
        box.Add(boxbuttons, 1, wx.ALIGN_CENTER, 3 )
        self.SetSizer(box)
        self.Centre()  
        
        #bindings
        self.Bind(wx.EVT_BUTTON, self.OnOk,self.ok,id=wx.ID_OK)
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)
        
        
    def detectCPU(self):
        """ detect the number of available CPU """
        if hasattr(os, "sysconf"):
            if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
                # Linux & Unix
                ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
                if isinstance(ncpus, int) and ncpus > 0:
                    return ncpus
            else:
                # OSX
                return int(os.popen2("sysctl -n hw.ncpu")[1].read())

        # Windows
        if os.environ.has_key("NUMBER_OF_PROCESSORS"):
            ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
            if ncpus > 0:
                return ncpus
        return 1
        
    def OnOk(self,e):
        pso_path=os.path.realpath(os.path.dirname(sys.argv[0]))
        print "> launching optimization with "+str(self.editcpu.GetValue())+" processors..."
        call='mpirun -n '+str(self.editcpu.GetValue())+' '+pso_path+'/POW.py DockSymmCircle pso.txt'
        subprocess.check_call(call,shell=True)
        os.remove('pso.txt')
        
        self.Close(True)

    def OnCancel(self,e):
        self.Close(True)

        
class IOPanel(wx.Panel):
    def __init__(self, parent):
        
        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)        
    
        self.dirname=''
        self.parent=parent
    
        #title       
        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="STRUCTURE")
        self.title.SetFont(font1)
        
        # neighbor type
        self.assemblystyle = ['rigid', 'flexible']
        self.lblassemblystyle = wx.StaticText(self, label="style:")
        self.editassemblystyle = wx.ComboBox(self, value=self.assemblystyle[0], size=(70, -1), choices=self.assemblystyle, style=wx.CB_DROPDOWN)
    
        #stoichiometry
        self.lbldegree=wx.StaticText(self, label="stoichiometry:")
        self.setdegree=wx.TextCtrl(self, value="",size=(35, -1))
    
        gs1 = wx.GridSizer(1, 4, 3, 3)
        gs1.AddMany([(self.lblassemblystyle, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editassemblystyle, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lbldegree, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.setdegree, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])
    
        #structure line
        self.lblstructure=wx.StaticText(self, label="structure:")
        self.setstructure=wx.TextCtrl(self, value="")
        self.structurebutton=wx.Button(self, label="select")
        
        #topology line
        self.lbltop=wx.StaticText(self, label="topology:")
        self.settop=wx.TextCtrl(self, value="")
        self.topbutton=wx.Button(self, label="select")
        self.settop.Enable(False)
        self.topbutton.Enable(False)

        #projection file
        self.lblproj=wx.StaticText(self, label="projection:")
        self.setproj=wx.TextCtrl(self, value="")
        self.projbutton=wx.Button(self, label="select")
        self.setproj.Enable(False)
        self.projbutton.Enable(False)

        #ratio
        self.lblratio=wx.StaticText(self, label="ratio:")
        self.setratio=wx.TextCtrl(self, value="0.9",size=(35, -1))
        self.setratio.Enable(False)
        
        #put structure, topology and projection stuff in a grid
        gs2 = wx.FlexGridSizer(3, 3, 3, 10)
        gs2.SetFlexibleDirection(wx.HORIZONTAL)
        gs2.AddGrowableCol(1,1)
        gs2.AddMany([(self.lblstructure, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setstructure, 1,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.structurebutton, 0, wx.ALIGN_LEFT, 3),
                    (self.lbltop, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3),
                    (self.settop, 1, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.topbutton, 0, wx.ALIGN_LEFT, 3),
                    (self.lblproj, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3),
                    (self.setproj, 1, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.projbutton, 0, wx.ALIGN_LEFT, 3),
                    (self.lblratio, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setratio, 0, wx.ALIGN_LEFT, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3)])

        self.Bind(wx.EVT_BUTTON, self.OnClickSelectStruct,self.structurebutton)
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectTop,self.topbutton)
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectProj,self.projbutton)    
        self.Bind(wx.EVT_COMBOBOX, self.OnStyleCombobox,self.editassemblystyle)   
        
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs1, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs2, 1, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()  

    def OnStyleCombobox(self,event):
        #if the style is rigid, disable everything, enable otherwise
        if self.editassemblystyle.Value=="rigid":
            self.activator=False
            #self.settop.Hide()
            #self.topbutton.Hide()
            #self.lbltop.Hide()
            #self.setproj.Hide()
            #self.projbutton.Hide()
            #self.lblproj.Hide()
            #self.setratio.Hide()        
            #self.lblratio.Hide()
        else:
            self.activator=True        
            #self.settop.Show()
            #self.topbutton.Show()
            #self.lbltop.Show()
            #self.setproj.Show()
            #self.projbutton.Show()
            #self.lblproj.Show()
            #self.setratio.Show()        
            #self.lblratio.Show()

        self.settop.Enable(self.activator)
        self.topbutton.Enable(self.activator)
        self.setproj.Enable(self.activator)
        self.projbutton.Enable(self.activator)
        self.setratio.Enable(self.activator)

    def OnClickSelectStruct(self,event):
        dlg = wx.FileDialog(self, "Choose your structure (pdb or trajectory)", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setstructure.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()

    def OnClickSelectTop(self,event):
        dlg = wx.FileDialog(self, "Choose your topology file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.settop.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()
    
    def OnClickSelectProj(self,event):
        dlg = wx.FileDialog(self, "Choose your projection file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setproj.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()
    
    
    def check(self):
        """ Verify parameters consistency"""
      
        #check stoichiometry value
        try:
            setdegree_f=int(self.setdegree.GetValue())
        except ValueError:
            self.parent.errorPopup("soichiometry should be an integer number greater than 1")
            return -1
        
        if setdegree_f<2:
            self.parent.errorPopup("assembly soichiometry should be greater or equal to 2")            
            return -1
      
        #check constraint file existence 
        fitfile=self.setstructure.GetValue()
        if fitfile == '':
            self.parent.errorPopup("structure not defined!")
            return -1 
        if os.path.isfile(fitfile)==False:
            self.parent.errorPopup("structure does not exist!")
            return -1 

        if self.editassemblystyle.GetValue()=="rigid":
            return 0
      
        #check topology file existence 
        fitfile=self.settop.GetValue()
        if fitfile == '':
            self.parent.errorPopup("when performing flexible assembly, a topology file should be provided (crd or dcd format)!")
            return -1 
        if os.path.isfile(fitfile)==False:
            self.parent.errorPopup("topology file does not exist!")
            return -1 

        #check projection file existence 
        projfile=self.setproj.GetValue()
        if projfile!='' and os.path.isfile(projfile)==False:
            self.parent.errorPopup("projections file does not exist!")
            return -1 

        if projfile=='':
            try:
                setratio_f=float(self.setratio.GetValue())
            except ValueError:
                self.parent.errorPopup("ratio value undefined! Expecting a number within 0 and 1.")
                return -1
            if setratio_f<0 or setratio_f>1:
                self.parent.errorPopup("ratio should be a number within 0 and 1!")            
                return -1

class ReceptorPanel(wx.Panel):
    def __init__(self, parent):

        self.dirname=''
        self.parent=parent

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)        
        
        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.t = wx.StaticText(self, label="RECEPTOR")
        self.t.SetFont(font1)
        # use repulsion field
        self.usereceptor = wx.CheckBox(self, label="")

        self.title = wx.FlexGridSizer(1, 2, 3, 5)
        self.title.SetFlexibleDirection(wx.HORIZONTAL)
        self.title.AddMany([(self.t, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT,3),
                    (self.usereceptor, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 3)])

        #structure line
        self.lblstructure=wx.StaticText(self, label="structure:")
        self.setstructure=wx.TextCtrl(self, value="")
        self.structurebutton=wx.Button(self, label="select")
        
        #box the group
        gs1 = wx.FlexGridSizer(1, 3, 3, 5)
        gs1.SetFlexibleDirection(wx.HORIZONTAL)
        gs1.AddGrowableCol(1,1)
        gs1.AddMany([(self.lblstructure, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setstructure, 1,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.structurebutton, 0, wx.ALIGN_LEFT, 3)])
        
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectReceptor,self.structurebutton)
        
        #top lables        
        self.lblx = wx.StaticText(self, label="rot. x")
        self.lbly = wx.StaticText(self, label="rot. y")
        self.lblz = wx.StaticText(self, label="rot. z")
        self.lblr = wx.StaticText(self, label="trans. z")

        #min boundaries
        self.lblmin = wx.StaticText(self, label="min :")
        self.minx = wx.TextCtrl(self, value="0", size=(40,-1))
        self.miny = wx.TextCtrl(self, value="0", size=(40,-1))
        self.minz = wx.TextCtrl(self, value="0", size=(40,-1))
        self.minr = wx.TextCtrl(self, value="", size=(40,-1))
        
        #max boundaries
        self.lblmax = wx.StaticText(self, label="max :")
        self.maxx = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxy = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxz = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxr = wx.TextCtrl(self, value="", size=(40,-1))
        
        #grid of PSO constants
        gs = wx.GridSizer(3, 5, 3, 3)
        gs.AddMany([(wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lblx, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lbly, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblz, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblr, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblmin, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minx, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.miny, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minz, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minr, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblmax, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxx, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxy, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxz, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxr, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        self.Bind(wx.EVT_CHECKBOX, self.OnUseReceptor,self.usereceptor)
        self.structurebutton.Enable(False)
        self.setstructure.Enable(False)
        self.minx.Enable(False)
        self.miny.Enable(False)
        self.minz.Enable(False)
        self.minr.Enable(False)
        self.maxx.Enable(False)
        self.maxy.Enable(False)
        self.maxz.Enable(False)
        self.maxr.Enable(False)


        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        #box.Add(self.usereceptor, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs1, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        #box.Add(self.repel, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()  
        
    def OnClickSelectReceptor(self,event):
        dlg = wx.FileDialog(self, "Choose your receptor file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setstructure.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


    def OnUseReceptor(self,event):
        #if the style is rigid, disable everything, enable otherwise
        if self.usereceptor.Value==False:
            self.activator=False
        else:
            self.activator=True        

        self.structurebutton.Enable(self.activator)
        self.setstructure.Enable(self.activator)
        self.minx.Enable(self.activator)
        self.miny.Enable(self.activator)
        self.minz.Enable(self.activator)
        self.minr.Enable(self.activator)
        self.maxx.Enable(self.activator)
        self.maxy.Enable(self.activator)
        self.maxz.Enable(self.activator)
        self.maxr.Enable(self.activator)


    def check(self):
        """ Verify parameters consistency"""
      
        if self.usereceptor.GetValue()==False:
            return 0
      
        #check constraint file existence 
        self.setstructure.GetValue()
        fitfile=self.setstructure.GetValue()
        if fitfile == '':
            self.parent.errorPopup("receptor structure not defined!")
            return -1 
        if os.path.isfile(fitfile)==False:
            self.parent.errorPopup("receptor structure does not exist!")
            return -1 

        #type consistency for min and max boundaries  
        try:
            minx_f=float(self.minx.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor x min boundary must be defined as a number within 0 and 360")
            return -1
   
        try:
            miny_f=float(self.miny.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor y min boundary must be defined as a number within 0 and 360")
            return -1
  
        try:
            minz_f=float(self.minz.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor z min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            minr_f=float(self.minr.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor radius min boundary must be defined as a real number")
            return -1
            
        try:
            maxx_f=float(self.maxx.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor x min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            maxy_f=float(self.maxy.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor y max boundary must be defined as a number within 0 and 360")
            return -1
  
        try:
            maxz_f=float(self.maxz.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor z max boundary must be defined as a number within 0 and 360")
            return -1
                 
        try:
            maxr_f=float(self.maxr.GetValue())
        except ValueError:
            self.parent.errorPopup("receptor radius max boundary must be defined as a real number")
            return -1    

        #value consistency for min and max boundaries
        if minx_f>maxx_f:
            self.parent.errorPopup("receptor in rotation around x axis, min boundary is greater than max!")
            return -1            

        if miny_f>maxy_f:
            self.parent.errorPopup("receptor in rotation around y axis, min boundary is greater than max!")
            return -1  
        
        if minz_f>maxz_f:
            self.parent.errorPopup("receptor in rotation around z axis, min boundary is greater than max!")
            return -1  

        if minr_f>maxr_f:
            self.parent.errorPopup("receptor in assembly radius definition, min boundary is greater than max!")
            return -1          


class ConstraintPanel(wx.Panel):
    def __init__(self, parent):

        self.dirname=''
        self.parent=parent

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)   
        
        self.font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)

        self.title = wx.StaticText(self, label="FITNESS")
        self.title.SetFont(self.font1)

        #fitness line
        self.lblconstraint=wx.StaticText(self, label="constraint:")
        self.setconstraint=wx.TextCtrl(self, value="")
        self.constraintbutton=wx.Button(self, label="select")
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectConstraint,self.constraintbutton)
        
        # boundary conditions
        self.lc = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
        self.lc.InsertColumn(0, 'ID')
        self.lc.InsertColumn(1, 'value')
        self.lc.InsertColumn(2, 'keyword')
        self.lc.SetColumnWidth(0, 50)

        #boundary conditions buttons box
        self.buttonsbox = wx.BoxSizer(wx.VERTICAL)
        self.buttonsbox.Add(wx.Button(self, 1, "Add"), 0 )
        self.buttonsbox.Add(wx.Button(self, 2, "Edit"), 0 )
        self.buttonsbox.Add(wx.Button(self, 3, "Remove"), 0 )

        #box the group
        gs = wx.FlexGridSizer(2, 3, 3, 5)
        #gs.SetFlexibleDirection(wx.HORIZONTAL)
        gs.AddGrowableCol(1,1)
        gs.AddGrowableRow(1,1)
        gs.AddMany([(self.lblconstraint, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3),
                    (self.setconstraint, 0, wx.ALIGN_TOP | wx.ALIGN_LEFT| wx.EXPAND, 3),
                    (self.constraintbutton, 0, wx.ALIGN_LEFT, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lc, 1, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.buttonsbox, 0, wx.ALIGN_TOP | wx.ALIGN_LEFT, 3)])

        # energy vs geometry mixing weight
        self.lblmix = wx.StaticText(self, label="mix coeff.:")
        self.editmix = wx.TextCtrl(self, value="0.2", size=(40,-1))

	#log filtering criteria
        self.lblfilter = wx.StaticText(self, label="filtering thresh.:")
        self.editfilter = wx.TextCtrl(self, value="0", size=(40,-1))

        gs2 = wx.GridSizer(1, 4, 3, 3)
        gs2.AddMany([(self.lblmix, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editmix, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblfilter, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editfilter, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        #draw title, boundaries and buttons list
        boundbox = wx.BoxSizer(wx.VERTICAL)
        boundbox.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        boundbox.Add(gs, 1, wx.EXPAND | wx.ALL, 3 )
        boundbox.Add(gs2, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(boundbox,1)
        self.SetSize((400,600))
        self.Centre()  

        #buttons bindings
        self.Bind(wx.EVT_BUTTON, self.OnBadd,id=1)
        self.Bind(wx.EVT_BUTTON, self.OnBedit,id=2)
        self.Bind(wx.EVT_BUTTON, self.OnBremove,id=3)


    def OnClickSelectConstraint(self,event):
        dlg = wx.FileDialog(self, "Choose your constraint file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setconstraint.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


    def OnBadd(self,event):

        dia=Boundaries(self)
        dia.ShowModal()
        dia.Destroy()

       
    def OnBedit(self,event):

        index = self.lc.GetFocusedItem()
        
        params=[]
        params.append(index)
        if index != -1 :
            params.append(self.lc.GetItem(index, 1).GetText())
            params.append(self.lc.GetItem(index, 2).GetText())
            params.append(self.lc.GetItem(index, 3).GetText())
            params.append(self.lc.GetItem(index, 4).GetText())
            params.append(self.lc.GetItem(index, 5).GetText())
        
            dia=Boundaries(self,params)
            dia.ShowModal()
            dia.Destroy()
                
    def OnBremove(self,event):
        index = self.lc.GetFocusedItem()
        self.lc.DeleteItem(index)
        
        #renumber items after deletion
        num_items = self.lc.GetItemCount()
        for x in xrange(index,num_items,1):
            self.lc.SetStringItem(x, 0, str(x+1))
           
    def check(self):
        #check constraint file existence 
        self.setconstraint.GetValue()
        fitfile=self.setconstraint.GetValue()
        if fitfile == '':
            self.parent.errorPopup("constraint file not defined!")
            return -1 
        if os.path.isfile(fitfile)==False:
            self.parent.errorPopup("constraint file does not exist!")
            return -1 

        l=self.lc.GetItemCount()
        if l == 0:
            self.parent.errorPopup("no constraints have been defined!")
            return -1        
        #test that the number of values outputted by the function and the number of lines is the same

        for x in xrange(0,l,1):      
            #check min and max
            m1=self.lc.GetItem(x, 1).GetText()
            if m1=='':
                self.parent.errorPopup(" constraint %s has no value set!"%(x+1))
                return -1 
            try:
                m1_f=float(m1)
            except ValueError:
                self.parent.errorPopup("constraint %s should be a number"%(x+1))
                return -1
            
        try:
            mix_f=float(self.editmix.GetValue())
        except ValueError:
            self.parent.errorPopup("mixing coefficient should be a number!")
            return -1           
        if mix_f<0 or mix_f>1:
            self.parent.errorPopup("mixing coefficient should be a number within 0 and 1")
            return -1


class BoundaryPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)        

        self.parent=parent
        
        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="BOUNDARIES")
        self.title.SetFont(font1)

        #top lables        
        self.lblx = wx.StaticText(self, label="rot. x")
        self.lbly = wx.StaticText(self, label="rot. y")
        self.lblz = wx.StaticText(self, label="rot. z")
        self.lblr = wx.StaticText(self, label="radius")

        #min boundaries
        self.lblmin = wx.StaticText(self, label="min :")
        self.minx = wx.TextCtrl(self, value="0", size=(40,-1))
        self.miny = wx.TextCtrl(self, value="0", size=(40,-1))
        self.minz = wx.TextCtrl(self, value="0", size=(40,-1))
        self.minr = wx.TextCtrl(self, value="", size=(40,-1))
        
        #max boundaries
        self.lblmax = wx.StaticText(self, label="max :")
        self.maxx = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxy = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxz = wx.TextCtrl(self, value="360", size=(40,-1))
        self.maxr = wx.TextCtrl(self, value="", size=(40,-1))
        
        #grid of PSO constants
        gs = wx.GridSizer(3, 5, 3, 3)
        gs.AddMany([(wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.lblx, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lbly, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblz, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblr, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblmin, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minx, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.miny, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minz, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.minr, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblmax, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxx, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxy, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxz, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.maxr, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        #box.Add(self.repel, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box,1)
        self.Centre()  

    def check(self):
        """ Verify parameters consistency"""
      
        #type consistency for min and max boundaries  
        try:
            minx_f=float(self.minx.GetValue())
        except ValueError:
            self.parent.errorPopup("x min boundary must be defined as a number within 0 and 360")
            return -1
   
        try:
            miny_f=float(self.miny.GetValue())
        except ValueError:
            self.parent.errorPopup("y min boundary must be defined as a number within 0 and 360")
            return -1
  
        try:
            minz_f=float(self.minz.GetValue() )
        except ValueError:
            self.parent.errorPopup("z min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            minr_f=float(self.minr.GetValue())
        except ValueError:
            self.parent.errorPopup("radius min boundary must be defined as a real number")
            return -1
            
        try:
            maxx_f=float(self.maxx.GetValue())
        except ValueError:
            self.parent.errorPopup("x min boundary must be defined as a number within 0 and 360")
            return -1

        try:
            maxy_f=float(self.maxy.GetValue())
        except ValueError:
            self.parent.errorPopup("y max boundary must be defined as a number within 0 and 360")
            return -1
  
        try:
            maxz_f=float(self.maxz.GetValue())
        except ValueError:
            self.parent.errorPopup("z max boundary must be defined as a number within 0 and 360")
            return -1
                 
        try:
            maxr_f=float(self.maxr.GetValue())
        except ValueError:
            self.parent.errorPopup("radius max boundary must be defined as a real number")
            return -1    

        #value consistency for min and max boundaries
        if minx_f>maxx_f:
            self.parent.errorPopup("in rotation around x axis, min boundary is greater than max!")
            return -1            

        if miny_f>maxy_f:
            self.parent.errorPopup("in rotation around y axis, min boundary is greater than max!")
            return -1  
        
        if minz_f>maxz_f:
            self.parent.errorPopup("in rotation around z axis, min boundary is greater than max!")
            return -1  

        if minr_f>maxr_f:
            self.parent.errorPopup("in assembly radius definition, min boundary is greater than max!")
            return -1          


     
class ProtocolPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER) 
        
        self.parent=parent
        
        #PSO behavior values     
        self.editneighsize=1
        self.editneigh="geographic"
        self.editinertiamin=0.3
        self.editinertiamax=0.7
        self.editpersonal=1.2
        self.editglobal=1.4
        
        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="PROTOCOL")
        self.title.SetFont(font1)

        # PSO steps
        self.lblsteps = wx.StaticText(self, label="Steps")
        self.editsteps = wx.TextCtrl(self, value="100", size=(40,-1))
        # particles
        self.lblpart = wx.StaticText(self, label="Particles")
        self.editpart = wx.TextCtrl(self, value="40", size=(40,-1))
        #protocol
        self.lblrep = wx.StaticText(self, label="Repetitions")
        self.editrep = wx.TextCtrl(self, value="1", size=(40,-1))

        self.detailsbutton=wx.Button(self, label="Details...")
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectDetails,self.detailsbutton)

        gs = wx.GridSizer(2, 4, 3, 3)
        gs.AddMany([(self.lblsteps, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblpart, 0,  wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblrep, 0,  wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (wx.Size(1,1), 0, wx.ALIGN_LEFT, 3),
                    (self.editsteps, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editpart, 0,  wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editrep, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.detailsbutton, 0, wx.ALIGN_CENTER | wx.ALIGN_CENTER_VERTICAL, 3)])
        
        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()         
        
    def OnClickSelectDetails(self,event):
        params=[self.editneigh,self.editneighsize,self.editinertiamin,self.editinertiamax,self.editpersonal,self.editglobal]
        dia=PSO(self,params)
        dia.ShowModal()
        dia.Destroy()
      
      
    def check(self):
        #check inertia
        imin=self.editinertiamin
        imax=self.editinertiamax
        try:
            imin_f=float(imin)
        except ValueError:
            self.parent.errorPopup("min inertia must be defined as a number within 0 and 1")
            return -1
        try:
            imax_f=float(imax)
        except ValueError:
            self.parent.errorPopup("max inertia must be defined as a number within 0 and 1")
            return -1           
        if imin_f<0:
            self.parent.errorPopup("min inertia should be greater than zero")
            return -1
        if imax_f>=1:
            self.parent.errorPopup("max inertia should be smaller than 1")
            return -1
        if imin_f>imax_f:
            self.parent.errorPopup("max inertia should be greater or equal than min inertia")
            return -1
        
        cn=self.editglobal
        try:
            cn_f=float(cn)
        except ValueError:
            self.parent.errorPopup("global weight should be a number")
            return -1
        if cn_f<0:
            self.parent.errorPopup("global weight should be positive")
            return -1
        
        cp=self.editpersonal
        try:
            cp_f=float(cp)
        except ValueError:
            self.parent.errorPopup("personal weight should be a number")
            return -1
        if cp_f<0:
            self.parent.errorPopup("personal weight should be positive")
            return -1
        
        psosteps=self.editsteps.GetValue()
        try:
            psosteps_f=int(psosteps)
        except ValueError:
            self.parent.errorPopup("number of PSO steps should be an integer")
            return -1
        if psosteps_f<0:
            self.parent.errorPopup("number of PSO steps should be greater than zero")
            return -1

        rep=self.editrep.GetValue()
        try:
            rep_f=int(rep)
        except ValueError:
            self.parent.errorPopup("number of PSO repetitions should be an integer")
            return -1
        if rep_f<0:
            self.parent.errorPopup("number of PSO repetitions should be greater than zero")
            return -1
        
        particles=self.editpart.GetValue()
        try:
            particles_f=int(particles)
        except ValueError:
            self.parent.errorPopup("number of particles should be an integer")
            return -1
        if particles_f<0:
            self.parent.errorPopup("number of particles should be greater than zero")
            return -1     

        nsize=self.editneighsize
        try:
            nsize_f=int(nsize)
        except ValueError:
            self.parent.errorPopup("neighborhood size should be an integer")
            return -1
        if nsize_f<0:
            self.parent.errorPopup("neighborhood size should be greater than zero")
            return -1
        if nsize_f>=particles_f:
            self.parent.errorPopup("neighborhood size should be smaller than number of particles")
            return -1
            
        ntype=self.editneigh
        if ntype!='geographic' and ntype!='indexed':
            self.parent.errorPopup("neighborhood type should be either geographic, either indexed")
            return -1 

        
        
class PSO(wx.Dialog):
    def __init__(self, parent,p):

        wx.Dialog.__init__(self, parent,title="PSO parametrization",size=(350,125))
        
        self.parent=parent
        
        # neighbor type
        self.neighList = ['geographic', 'indexed']
        self.lblneigh = wx.StaticText(self, label="Neighbor style : ")
        self.editneigh = wx.ComboBox(self, value=str(p[0]), size=(80, -1), choices=self.neighList, style=wx.CB_DROPDOWN)
        # neighborhood size
        self.lblneighsize = wx.StaticText(self, label="Neighbors nb. :")
        self.editneighsize = wx.TextCtrl(self, value=str(p[1]), size=(40,-1))

        self.lblinertiamin = wx.StaticText(self, label="inertia min :")
        self.editinertiamin = wx.TextCtrl(self, value=str(p[2]), size=(40,-1))
        self.lblinertiamax = wx.StaticText(self, label="inertia max :")
        self.editinertiamax = wx.TextCtrl(self, value=str(p[3]), size=(40,-1))

        self.lblpersonal = wx.StaticText(self, label="personal weight :")
        self.editpersonal = wx.TextCtrl(self, value=str(p[4]), size=(40,-1))
        self.lblglobal = wx.StaticText(self, label="global weight :")
        self.editglobal = wx.TextCtrl(self, value=str(p[5]), size=(40,-1))

        # use repulsion field
        #self.repel = wx.CheckBox(self, label="use repulsion field?")
        
        #grid of PSO constants
        gs = wx.GridSizer(3, 4, 3, 3)
        gs.AddMany([(self.lblneigh, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editneigh, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblneighsize, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editneighsize, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblinertiamin, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editinertiamin, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblinertiamax, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editinertiamax, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblpersonal, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editpersonal, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblglobal, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editglobal, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        
        self.okbutton=wx.Button(self, label="OK!")
        self.Bind(wx.EVT_BUTTON, self.OnOk,self.okbutton,id=wx.ID_OK)
        self.cancel =wx.Button(self, label="Cancel")
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)
        self.buttonbox = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonbox.Add(self.okbutton, 0, wx.RIGHT, 3 )
        self.buttonbox.Add(self.cancel, 0, wx.LEFT, 3 )

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.buttonbox, 1, wx.CENTER, 3 )
        self.SetSizer(box)
        self.Centre()  

    def OnOk(self,e):
        self.parent.editneighsize=self.editneighsize.GetValue()
        self.parent.editneigh=self.editneigh.GetValue()
        self.parent.editinertiamin=self.editinertiamin.GetValue()
        self.parent.editinertiamax=self.editinertiamax.GetValue()
        self.parent.editpersonal=self.editpersonal.GetValue()
        self.parent.editglobal=self.editglobal.GetValue()
        
        self.Close(True)  # Close the frame.

    def OnCancel(self,e):
        self.Close(True)


class Boundaries(wx.Dialog):
    def __init__(self, parent,params=-1):
        wx.Dialog.__init__(self, parent,title="Boundaries Editor",size=(300,75))
        
        self.params=params
        self.parent=parent
        
        self.lblval = wx.StaticText(self, label="value:")
        self.editval = wx.TextCtrl(self, value="", size=(40,-1))
        self.lblkeyword = wx.StaticText(self, label="keyword:")
        self.editkeyword = wx.TextCtrl(self, value="", size=(75,-1))

        if params!=-1:
            self.editval.SetValue(self.params[1])
            self.editkeyword.SetValue(self.params[2])

        gs = wx.GridSizer(1, 4, 3, 3)
        gs.AddMany([(self.lblval, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editval, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblkeyword, 1,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editkeyword, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        self.okbutton=wx.Button(self, label="OK!")
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.okbutton,id=wx.ID_OK)
        self.cancel =wx.Button(self, label="Cancel")
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)
        self.buttonbox = wx.BoxSizer(wx.HORIZONTAL)
        self.buttonbox.Add(self.okbutton, 0, wx.RIGHT, 3 )
        self.buttonbox.Add(self.cancel, 0, wx.LEFT, 3 )

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(gs, 0, wx.EXPAND, 3 )
        box.Add(self.buttonbox, 1, wx.CENTER, 3 )
        self.SetSizer(box)
        self.Centre()

    def OnOk(self,e):
        
        if self.params!=-1:
            num_items=self.params[0]            
        else:
            num_items = self.parent.lc.GetItemCount()
            self.parent.lc.InsertStringItem(num_items, str(num_items+1))
            
        self.parent.lc.SetStringItem(num_items, 1, self.editval.GetValue())
        self.parent.lc.SetStringItem(num_items, 2, self.editkeyword.GetValue())
        
        self.Close(True)  # Close the frame.

    def OnCancel(self,e):
        self.Close(True)  # Close the frame.


app = wx.App(False)
frame = MainWindow(None, "POW Assembly of Symmetrical Assemblies")
frame.Show()
app.MainLoop()
