#! /usr/bin/env python

import os, sys
import wx
import subprocess


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

        # creating the box sizer and addint the buttons to it!
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        buttons.Add(self.verify, 0 )
        buttons.Add(self.launch, 0 )

        #panels
        self.I=IOPanel(self) # those are the panels as form of classes
        self.B=BoundaryPanel(self)
        self.P=PSOPanel(self)
        self.S=ProtocolPanel(self)

        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.I, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.B, 1, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.P, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.S, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(buttons, 0, wx.ALIGN_CENTER)
        self.SetSizer(box)
        self.SetSize((450,600))
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
        dlg = wx.MessageDialog(self, " pso4lbm v 2.0 \n Matteo Degiacomi \n 2011", "About pso4lbm", wx.OK)
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
        """ Launch call to PDF reference manual"""
        pass

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

    def nullify(self):
        """ Reset all the form fields"""
        self.S.editsteps.SetValue('')
        self.S.editpart.SetValue('')
        self.S.editrep.SetValue('')
        self.P.editneighsize.SetValue('')
        self.P.editneigh.SetValue('')
        self.P.editinertiamax.SetValue('')
        self.P.editinertiamin.SetValue('')
        self.P.editpersonal.SetValue('')
        self.P.editglobal.SetValue('')
        self.I.setlog.SetValue('')
        self.I.setfitness.SetValue('')
        self.P.repel.SetValue(False)

        for __ in xrange(0,self.B.lc.GetItemCount(),1):
            self.B.lc.DeleteItem(0)


    def parse(self,infile):
        """ Input file parser"""

        self.nullify()

        f=open(infile, 'r')
        line = f.readline()

        low=''
        high=''
        max_vel=''
        boundary_type=''
        keyword=''
        while line:
            w = line.split()
            if len(w) > 0 :
                if str(w[0])=='repeat' : self.S.editrep.SetValue(str(w[1]))
                if str(w[0])=='steps' : self.S.editsteps.SetValue(str(w[1]))
                if str(w[0])=='particles' : self.S.editpart.SetValue(str(w[1]))
                if str(w[0])=='neighborSize' : self.P.editneighsize.SetValue(str(w[1]))
                if str(w[0])=='neighborType' : self.P.editneigh.SetValue(str(w[1]))
                if str(w[0])=='inertiaMax' : self.P.editinertiamax.SetValue(str(w[1]))
                if str(w[0])=='inertiaMin' : self.P.editinertiamin.SetValue(str(w[1]))
                if str(w[0])=='cp' : self.P.editpersonal.SetValue(str(w[1]))
                if str(w[0])=='cn' : self.P.editglobal.SetValue(str(w[1]))
                if str(w[0])=='output' : self.I.setlog.SetValue(str(w[1]))
                if str(w[0])=='fitnessFile' : self.I.setfitness.SetValue(str(w[1]))
                if str(w[0])=='repulsion' :
                    if str(w[1]) == 'on':
                        self.P.repel.SetValue(True)
                    else:
                        self.P.repel.SetValue(False)
                if str(w[0])=='boundaryMin' : low=w[1:len(w)]
                if str(w[0])=='boundaryMax' : high=w[1:len(w)]
                if str(w[0])=='velocityMax' : max_vel=w[1:len(w)]
                if str(w[0])=='boundaryType' : boundary_type=w[1:len(w)]
                if str(w[0])=='keyword' : keyword=w[1:len(w)]
            line = f.readline()
        f.close()
        #insert data from boundary conditions list
        for x in xrange(0,len(low),1):
            self.B.lc.InsertStringItem(x, str(x+1))
            if low!='' and str(low[x]) != "NaN":
                self.B.lc.SetStringItem(x, 1, str(low[x]))
            if high!='' and str(high[x]) != "NaN":
                self.B.lc.SetStringItem(x, 2, str(high[x]))
            if max_vel!='' and str(max_vel[x]) != "NaN":
                self.B.lc.SetStringItem(x, 3, str(max_vel[x]))
            if boundary_type!='' and str(boundary_type[x]) != "NaN":
                b=''
                if str(boundary_type[x]) == 'p':
                    b='periodic'
                elif str(boundary_type[x]) == 'r':
                    b='repulsive'
                self.B.lc.SetStringItem(x, 4, b)
            if keyword!='' and str(keyword[x]) != "NaN":
                self.B.lc.SetStringItem(x, 5, str(keyword[x]))

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

        r=self.check()
        if r!=-1:
            self.errorPopup("Setup OK!")

    def check(self):
        """ Verify parameters consistency"""

        fitfile=self.I.setfitness.GetValue()
        if fitfile == '':
            self.errorPopup("fitness file not defined!")
            return -1
        if os.path.isfile(fitfile)==False:
            self.errorPopup("fitness file does not exist!")
            return -1

        if self.I.setlog.GetValue()== '':
            self.errorPopup("log file not defined!")
            return -1

        l=self.B.lc.GetItemCount()
        if l == 0:
            self.errorPopup("no boundaries have been defined!")
            return -1

        for x in xrange(0,l,1):
            #check min and max
            m1=self.B.lc.GetItem(x, 1).GetText()
            m2=self.B.lc.GetItem(x, 2).GetText()
            if m1=='':
                self.errorPopup("boundary %s has no min set!"%(x+1))
                return -1
            if m2=='':
                self.errorPopup("boundary %s has no max set!"%(x+1))
                return -1
            try:
                m1_f=float(m1)
            except ValueError:
                self.errorPopup("min in boundary %s should be a number"%(x+1))
                return -1
            try:
                m2_f=float(m2)
            except ValueError:
                self.errorPopup("max in boundary %s should be a number"%(x+1))
                return -1
            if m2_f<m1_f:
                self.errorPopup("boundary %s inconsistent: min is greater than max!"%(x+1))
                return -1

            #check velocity
            v=self.B.lc.GetItem(x, 3).GetText()
            if v=='':
                v=0
            try:
                v_f=float(v)
            except ValueError:
                self.errorPopup("velocity in boundary %s should be a number"%(x+1))
                return -1
            if v_f<0:
                self.errorPopup("velocity in boundary %s should be positive!"%(x+1))
                return -1
            if v_f>m2_f-m1_f:
                self.errorPopup("velocity in boundary %s is too high! Max admitted velocity is %s"%(x+1,float(m2)-float(m1)))
                return -1

            #check neighborhood
            t=self.B.lc.GetItem(x, 4).GetText()
            if t!='repulsive' and t!='periodic' and t!='':
                self.errorPopup("type in boundary %s unknown"%(x+1))
                return -1

        #check inertia
        imin=self.P.editinertiamin.GetValue()
        imax=self.P.editinertiamax.GetValue()
        try:
            imin_f=float(imin)
        except ValueError:
            self.errorPopup("min inertia must be defined as a number within 0 and 1")
            return -1
        try:
            imax_f=float(imax)
        except ValueError:
            self.errorPopup("max inertia must be defined as a number within 0 and 1")
            return -1
        if imin_f<0:
            self.errorPopup("min inertia should be greater than zero")
            return -1
        if imax_f>=1:
            self.errorPopup("max inertia should be smaller than 1")
            return -1
        if imin_f>imax_f:
            self.errorPopup("max inertia should be greater or equal than min inertia")
            return -1

        cn=self.P.editglobal.GetValue()
        try:
            cn_f=float(cn)
        except ValueError:
            self.errorPopup("global weight should be a number")
            return -1
        if cn_f<0:
            self.errorPopup("global weight should be positive")
            return -1

        cp=self.P.editpersonal.GetValue()
        try:
            cp_f=float(cp)
        except ValueError:
            self.errorPopup("personal weight should be a number")
            return -1
        if cp_f<0:
            self.errorPopup("personal weight should be positive")
            return -1

        psosteps=self.S.editsteps.GetValue()
        try:
            psosteps_f=int(psosteps)
        except ValueError:
            self.errorPopup("number of PSO steps should be an integer")
            return -1
        if psosteps_f<0:
            self.errorPopup("number of PSO steps should be greater than zero")
            return -1

        rep=self.S.editrep.GetValue()
        try:
            rep_f=int(rep)
        except ValueError:
            self.errorPopup("number of PSO repetitions should be an integer")
            return -1
        if rep_f<0:
            self.errorPopup("number of PSO repetitions should be greater than zero")
            return -1

        particles=self.S.editpart.GetValue()
        try:
            particles_f=int(particles)
        except ValueError:
            self.errorPopup("number of particles should be an integer")
            return -1
        if particles_f<0:
            self.errorPopup("number of particles should be greater than zero")
            return -1

        nsize=self.P.editneighsize.GetValue()
        try:
            nsize_f=int(nsize)
        except ValueError:
            self.errorPopup("neighborhood size should be an integer")
            return -1
        if nsize_f<0:
            self.errorPopup("neighborhood size should be greater than zero")
            return -1
        if nsize_f>=particles_f:
            self.errorPopup("neighborhood size should be smaller than number of particles")
            return -1

        ntype=self.P.editneigh.GetValue()
        if ntype!='geographic' and ntype!='indexed':
            self.errorPopup("neighborhood type should be either geographic, either indexed")
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

        ffile=str(self.I.setfitness.GetValue())
        if ffile != '':
            f.write("fitnessFile "+ffile+"\n")

        lfile=str(self.I.setlog.GetValue())
        if lfile != '':
            f.write("output "+lfile+"\n")

        min=''
        max=''
        vel=''
        type=''
        keyword=''
        for x in xrange(0,self.B.lc.GetItemCount(),1):
            #min and max
            m1=self.B.lc.GetItem(x, 1).GetText()
            if m1=='':
                m1="NaN"
            min+=str(m1)+" "
            m2=self.B.lc.GetItem(x, 2).GetText()
            if m2=='':
                m2="NaN"
            max+=str(m2)+" "
            #velocity
            v=self.B.lc.GetItem(x, 3).GetText()
            if v =='':
                if estimate == False:
                    v="NaN"
                else:
                    v=float(m2)-float(m1)
            vel+=str(v)+" "
            #type
            if self.B.lc.GetItem(x, 4).GetText() == 'repulsive':
                type+=str('r')+" "
            else:
                type+=str('p')+" "
            #keyword
            k=self.B.lc.GetItem(x, 5).GetText()
            if k =='':
                k='NaN'
            else:
                keyword+=str(k)+" "

        if self.P.repel.GetValue() == True:
            r='on'
        else:
            r='off'

        if str(min)!='':
            f.write("boundaryMin "+str(min)+"\n")
        if str(max)!='':
            f.write("boundaryMax "+str(max)+"\n")
        if str(vel)!='':
            f.write("velocityMax "+str(vel)+"\n")
        if str(type)!='':
            f.write("boundaryType "+str(type)+"\n")
        if str(keyword)!='':
            f.write("keyword "+str(keyword)+"\n")
        if str(self.P.editinertiamin.GetValue())!='':
            f.write("inertiaMin "+self.P.editinertiamin.GetValue()+"\n")
        if str(self.P.editinertiamax.GetValue())!='':
            f.write("inertiaMax "+self.P.editinertiamax.GetValue()+"\n")
        if str(self.P.editglobal.GetValue())!='':
            f.write("cn "+self.P.editglobal.GetValue()+"\n")
        if str(self.P.editpersonal.GetValue())!='':
            f.write("cp "+self.P.editpersonal.GetValue()+"\n")
        if str(self.S.editrep.GetValue())!='':
            f.write("repeat "+self.S.editrep.GetValue()+"\n")
        if str(self.S.editsteps.GetValue())!='':
            f.write("steps "+self.S.editsteps.GetValue()+"\n")
        if str(self.S.editpart.GetValue())!='':
            f.write("particles "+self.S.editpart.GetValue()+"\n")
        if str(self.P.editneighsize.GetValue())!='':
            f.write("neighborSize "+self.P.editneighsize.GetValue()+"\n")
        if str(self.P.editneigh.GetValue())!='':
            f.write("neighborType "+self.P.editneigh.GetValue()+"\n")
        f.write("repulsion "+r+"\n")

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
        call='mpiexec -np '+str(self.editcpu.GetValue())+' '+pso_path+'/POW.py Function pso.txt'
        subprocess.check_call(call,shell=True)
        os.remove('pso.txt')

        self.Close(True)

    def OnCancel(self,e):
        self.Close(True)


class IOPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        self.dirname=''

        #title
        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="I/O FILES")
        self.title.SetFont(font1)

        #fitness line
        self.lblfitness=wx.StaticText(self, label="Fitness file:")
        self.setfitness=wx.TextCtrl(self, value="")
        self.fitbutton=wx.Button(self, label="select")

        #log line
        self.lbllog=wx.StaticText(self, label="Log file:")
        self.setlog=wx.TextCtrl(self, value="log.txt")
        self.logbutton=wx.Button(self, label="select")

        #put log and fitness stuff in a grid
        gs = wx.FlexGridSizer(2, 3, 20, 20) # flexGridSzer (self, rows, cols, vgap, hgap) info at http://wxpython.org/docs/api/wx.FlexGridSizer-class.html
        gs.SetFlexibleDirection(wx.HORIZONTAL) # specifies that rows are resized dynamically
        gs.AddGrowableCol(1,1) # specifies that the column 1(starting from 0) can be regrown, ids specified below!
        gs.AddMany([(self.lblfitness, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT,3), # here insert them in order plaase
                    (self.setfitness, 1,wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.fitbutton, 0, wx.ALIGN_LEFT, 3),
                    (self.lbllog, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 3),
                    (self.setlog, 1, wx.ALIGN_LEFT | wx.EXPAND, 3),
                    (self.logbutton, 0, wx.ALIGN_LEFT, 3)])

        self.Bind(wx.EVT_BUTTON, self.OnClickSelectFit,self.fitbutton)
        self.Bind(wx.EVT_BUTTON, self.OnClickSelectLog,self.logbutton)

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL) # more info about box sizing at http://zetcode.com/wxpython/layout/
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 1, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()


    def OnClickSelectFit(self,event):
        dlg = wx.FileDialog(self, "Choose your fitness function file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setfitness.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()

    def OnClickSelectLog(self,event):
        dlg = wx.FileDialog(self, "Choose your log file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.setlog.SetValue(os.path.join(self.dirname, self.filename))
        dlg.Destroy()


class BoundaryPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        self.font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)

        self.title = wx.StaticText(self, label="BOUNDARY CONDITIONS")
        self.title.SetFont(self.font1)

        # boundary conditions
        self.lc = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
        self.lc.InsertColumn(0, 'ID')
        self.lc.InsertColumn(1, 'min pos')
        self.lc.InsertColumn(2, 'max pos')
        self.lc.InsertColumn(3, 'max vel')
        self.lc.InsertColumn(4, 'type')
        self.lc.InsertColumn(5, 'keyword')
        self.lc.SetColumnWidth(0, 35)

        #boundary conditions buttons box
        buttonsbox = wx.BoxSizer(wx.HORIZONTAL)
        buttonsbox.Add(wx.Button(self, 1, "Add"), 0 )
        buttonsbox.Add(wx.Button(self, 2, "Edit"), 0 )
        buttonsbox.Add(wx.Button(self, 3, "Clone"), 0 )
        buttonsbox.Add(wx.Button(self, 4, "Remove"), 0 )

        #draw title, boundaries and buttons list
        boundbox = wx.BoxSizer(wx.VERTICAL)
        boundbox.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        boundbox.Add(self.lc, 1, wx.EXPAND | wx.ALL, 3 )
        boundbox.Add(buttonsbox, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(boundbox,1)
        self.SetSize((400,600))
        self.Centre()

        #buttons bindings
        self.Bind(wx.EVT_BUTTON, self.OnBadd,id=1)
        self.Bind(wx.EVT_BUTTON, self.OnBedit,id=2)
        self.Bind(wx.EVT_BUTTON, self.OnClone,id=3)
        self.Bind(wx.EVT_BUTTON, self.OnBremove,id=4)


    # creatin the boundary box from the class defined below
    def OnBadd(self,event):

        dia=Boundaries(self)
        dia.ShowModal() # making it principal on the screen
        dia.Destroy() # destroy it when done


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

    def OnClone(self,event):
        num_items = self.lc.GetItemCount()
        index = self.lc.GetFocusedItem()

        if index != -1 :
            min=self.lc.GetItem(index, 1).GetText()
            max=self.lc.GetItem(index, 2).GetText()
            vel=self.lc.GetItem(index, 3).GetText()
            type=self.lc.GetItem(index, 4).GetText()
            keyword=self.lc.GetItem(index, 5).GetText()

            self.lc.InsertStringItem(num_items, str(num_items+1))
            self.lc.SetStringItem(num_items, 1, min)
            self.lc.SetStringItem(num_items, 2, max)
            self.lc.SetStringItem(num_items, 3, vel)
            self.lc.SetStringItem(num_items, 4, type)
            self.lc.SetStringItem(num_items, 5, keyword)

    def OnBremove(self,event):
        index = self.lc.GetFocusedItem()
        self.lc.DeleteItem(index)

        #renumber items after deletion
        num_items = self.lc.GetItemCount()
        for x in xrange(index,num_items,1):
            self.lc.SetStringItem(x, 0, str(x+1))


class PSOPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="PSO BEHAVIOR")
        self.title.SetFont(font1)

        # neighbor type
        self.neighList = ['geographic', 'indexed']
        self.lblneigh = wx.StaticText(self, label="Neighbor style : ")
        self.editneigh = wx.ComboBox(self, value=self.neighList[0], size=(110, -1), choices=self.neighList, style=wx.CB_DROPDOWN)
        # neighborhood size
        self.lblneighsize = wx.StaticText(self, label="Neighbors nb. :")
        self.editneighsize = wx.TextCtrl(self, value="1", size=(40,-1))

        self.lblinertiamin = wx.StaticText(self, label="inertia min :")
        self.editinertiamin = wx.TextCtrl(self, value="0.3", size=(40,-1))
        self.lblinertiamax = wx.StaticText(self, label="inertia max :")
        self.editinertiamax = wx.TextCtrl(self, value="0.7", size=(40,-1))

        self.lblpersonal = wx.StaticText(self, label="personal weight :")
        self.editpersonal = wx.TextCtrl(self, value="1.2", size=(40,-1))
        self.lblglobal = wx.StaticText(self, label="global weight :")
        self.editglobal = wx.TextCtrl(self, value="1.4", size=(40,-1))

        # use repulsion field
        self.repel = wx.CheckBox(self, label="use repulsion field?")

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

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(self.repel, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()

class ProtocolPanel(wx.Panel):
    def __init__(self, parent):

        wx.Panel.__init__(self, parent,style=wx.SIMPLE_BORDER)

        font1 = wx.Font(10, wx.NORMAL, wx.ITALIC, wx.BOLD)
        self.title = wx.StaticText(self, label="PROTOCOL")
        self.title.SetFont(font1)

        # PSO steps
        self.lblsteps = wx.StaticText(self, label="Steps :")
        self.editsteps = wx.TextCtrl(self, value="100", size=(40,-1))
        # particles
        self.lblpart = wx.StaticText(self, label="Particles :")
        self.editpart = wx.TextCtrl(self, value="40", size=(40,-1))
        #protocol
        self.lblrep = wx.StaticText(self, label="Repetitions :")
        self.editrep = wx.TextCtrl(self, value="1", size=(40,-1))

        gs = wx.GridSizer(2, 4, 3, 3)
        gs.AddMany([(self.lblsteps, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editsteps, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblpart, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editpart, 0,  wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.lblrep, 0,  wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, 3),
                    (self.editrep, 0, wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL, 3)])

        #draw title, and grid in a vertical box
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.title, 0, wx.EXPAND | wx.ALL, 3 )
        box.Add(gs, 0, wx.EXPAND | wx.ALL, 3 )
        self.SetSizer(box)
        self.Centre()


class Boundaries(wx.Dialog):
    def __init__(self, parent,params=-1):
        wx.Dialog.__init__(self, parent,title="Boundaries Editor",size=(300,160))

        self.params=params
        self.parent=parent

        self.y=10
        self.lblmin = wx.StaticText(self, label="min:", pos=(20,self.y+3))
        self.editmin = wx.TextCtrl(self, value="", pos=(70,self.y), size=(40,-1))
        self.lblmax = wx.StaticText(self, label="max :", pos=(165,self.y+3))
        self.editmax = wx.TextCtrl(self, value="", pos=(200,self.y), size=(40,-1))

        self.y+=30

        # max vel
        self.lblnvel = wx.StaticText(self, label="max vel.:", pos=(20,self.y+3))
        self.editvel = wx.TextCtrl(self, value="", pos=(70,self.y), size=(40,-1))
        #boundary type
        self.typeList = ['periodic', 'repulsive']
        self.lbltype = wx.StaticText(self, label="boundary type : ", pos=(115, self.y+3))
        self.editneigh = wx.ComboBox(self, value=self.typeList[0], pos=(200, self.y), size=(80, -1), choices=self.typeList, style=wx.CB_DROPDOWN)

        self.y+=30

        self.lblkeyword = wx.StaticText(self, label="keyword:", pos=(20,self.y+3))
        self.editkeyword = wx.TextCtrl(self, value="", pos=(80, self.y), size=(200,-1))

        self.y+=30

        self.add =wx.Button(self, label="OK!", pos=(20, self.y))
        self.Bind(wx.EVT_BUTTON, self.OnOk,self.add,id=wx.ID_OK)
        self.cancel =wx.Button(self, label="Cancel", pos=(90, self.y))
        self.Bind(wx.EVT_BUTTON, self.OnCancel,self.cancel,id=wx.ID_CANCEL)

        if params!=-1:
            self.editmin.SetValue(self.params[1])
            self.editmax.SetValue(self.params[2])
            self.editvel.SetValue(self.params[3])
            self.editneigh.SetValue(self.params[4])
            self.editkeyword.SetValue(self.params[5])


    def OnOk(self,e):

        if self.params!=-1:
            num_items=self.params[0]
        else:
            num_items = self.parent.lc.GetItemCount()
            self.parent.lc.InsertStringItem(num_items, str(num_items+1))

        self.parent.lc.SetStringItem(num_items, 1, self.editmin.GetValue())
        self.parent.lc.SetStringItem(num_items, 2, self.editmax.GetValue())
        self.parent.lc.SetStringItem(num_items, 3, self.editvel.GetValue())
        self.parent.lc.SetStringItem(num_items, 4, self.editneigh.GetValue())
        self.parent.lc.SetStringItem(num_items, 5, self.editkeyword.GetValue())

        self.Close(True)  # Close the frame.

    def OnCancel(self,e):
        self.Close(True)  # Close the frame.


app = wx.App(False)
frame = MainWindow(None, "pso4lbm")
frame.Show()
app.MainLoop()
