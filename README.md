# Dislocation-Density
Crystal dislocation code

/*
    -- Main reference implementation "Nonlocal elasticity tensors in dislocation and disclination cores"
	  
    Authors to main reference: 
            Capolungo, Laurent
		        Gbemou, kodjovi
		        Taupin, Vincent
		        Fressengeas, Claude
    Info:        
		      Provided by the author(s) and the Los Alamos National Laboratory (2017-04-14).
		      To be published in: Journal of the Mechanics and Physics of Solids
		      DOI to publisher's version: 10.1016/j.jmps.2017.01.003
		      Permalink to record: http://permalink.lanl.gov/object/view?what=info:lanl-repo/lareport/LA-UR-17-20035

	  --> terminal compilation with g++ (--> = enter key):
		--> cd <current directory> --> g++ dislocationMain.cpp tensorMaps.cpp --> ./a.exe
    -- interfaces have extension (.h) and Implementations hav extension (.cpp).
    -- see data for data[i][j] for parameter decloration and instalization
    -- alpha must not be zero as it denotes the curl of the plastic destortion tensor U (See Eq.(15))
    -- T_i must be a unit vector as well as e_1,e_2, and e_3 --> add check method before runnning algo
*/

/*
 Author of algorithm:       Roy Burson
              
///
*/
