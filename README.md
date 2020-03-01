<!-- It seems that github simply ignores the "style" tags within div tags... so try something different -->
<p align=center>
<a href="http://icarus.lngs.infn.it"><img src="http://icarus.lngs.infn.it/img/n3.jpg" alt="ICARUS Experiment" style="border:0"></a>
</p>

<h1 align=center><font color="blue"><font size="7">ICARUS Signal Processing</font></font></h1><br>
<p align=center>
<font color="gray"><font size="3">A repository aimed at collecting common tools for handling waveforms/ROI/hit finding <br>with the ICARUS LAr TPC</font></font><br>
</p>


<h2><font color="blue"><font size="5">General Directory Structure</font></font></h2>
<ul>
    <li><b>icarus_signal_processing</b> - This folder contains our waveform processing tools</li>
    <li><b>Makefile</b> - This folder contains make file definitions for building without mrb</li>
    <li><b>test</b> - TODO - create some test programs and place in this folder</li>
    <li><b>ups</b> - This folder contains the product definitions for ups</li>
</ul>

<h2><font color="blue"><font size="5">Installing</font></font></h2>
<ul>
    <li>Create a working area on your machine where you will copy the input files and set up this repository</li>
    <li>Navigate to that folder
    <li>git clone https://github.com/SFBayLaser/icarus_signal_processing </li>
    <li>Two modes:</li>
        <ul>
            <li>Follow the standard procedure for building LArSoft packages using mrb</li>
            <li>Use make to build with the included GNuMakefile</li>
            <ul>
                <li>Note that you need to make sure you have ROOTSYS defined</li>
            </ul>
        </ul>
</ul>

<h2><font color="blue"><font size="5">Problems/Complaints/Hate mail</font></font></h2>
<ul>
    <li>Currently, all questions and/or complaints should be sent to <a href="mailto:usher@slac.stanford.edu">Tracy Usher</a></li>
</ul>

