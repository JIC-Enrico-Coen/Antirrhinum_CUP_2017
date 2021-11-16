function m = gpt_revisedcup_20170503( m )
%m = gpt_revisedcup_20170503( m )
%   Morphogen interaction function.
%   Written at 2017-07-13 16:48:06.
%   GFtbox revision 5501, 2017-06-06 09:41.

% The user may edit any part of this function between delimiters
% of the form "USER CODE..." and "END OF USER CODE...".  The
% delimiters themselves must not be moved, edited, deleted, or added.

    if isempty(m), return; end

    fprintf( 1, '%s found in %s\n', mfilename(), which(mfilename()) );

    try
        m = local_setproperties( m );
    catch
    end

    setGlobals();
    realtime = m.globalDynamicProps.currenttime;
    dt = m.globalProps.timestep;

%%% USER CODE: INITIALISATION

if (Steps(m)==0) && m.globalDynamicProps.doinit % First iteration
    
    
    m = leaf_setproperty( m, 'rectifyverticals', 1 );
    
    % Set up names for variant models.  Useful for running multiple models on a cluster.
    m.userdata.ranges.modelname.range {1} = 'isotropic growth Figure 4B';
    m.userdata.ranges.modelname.range {2} = 'isotropic growth inhibited at junction Figure 4C';
    m.userdata.ranges.modelname.range {3} = 'isotropic growth inhibited at sinus Figure 4D';
    m.userdata.ranges.modelname.range {4} = 'isotropic growth inhibited at sinus and base Figure 4E';
    m.userdata.ranges.modelname.range {6} = 'anisotropic growth Figure 4F';
    m.userdata.ranges.modelname.range {7} = 'anisotropic growth inhibited at junction Figure 4G';
    m.userdata.ranges.modelname.range {8} = 'anisotropic growth inhibited at sinus Figure 4H';
    m.userdata.ranges.modelname.range {9} = 'anisotropic growth inhibited at sinus and base Figure 4I';
    m.userdata.ranges.modelname.range {10} = 'growth promoted by med Figure 4K';
    m.userdata.ranges.modelname.range {11} = 'growth promoted by med and inhibited at junction Figure 4L';
    m.userdata.ranges.modelname.range {12} = 'growth promoted by med and inhibited at sinus Figure 4M';
    m.userdata.ranges.modelname.range {13} = 'growth promoted by med and selectively inhibited at sinus Figure 4N';
    m.userdata.ranges.modelname.index = 11;
    
    
    m.globalProps.timestep=2.5; %each step relates to 2.5 hours in real time
    
    % This is the code for creating a grid of cells around the cylinder.
    % However, it is not necessary to call this on every run, as the
    % initial mesh has been saved with the cells created.
    if ~hasNonemptySecondLayer( m )
        m = leaf_makesecondlayer( m, 'mode', 'cylindergrid', 'relarea', 0.0025 );  % This gives cells stacked 4 high.
        m = leaf_setproperty( m, 'bioAalpha', 0 );  % Make the cell grid transparent.
    end
end
modelname = m.userdata.ranges.modelname.range{m.userdata.ranges.modelname.index};
m = leaf_plotoptions( m, 'bioAalpha', 0 );

m.globalProps.timestep=2.5;% CLUSTER

disp(sprintf('\nRunning %s model %s\n',mfilename, modelname));


% set colour of polariser gradient arrows
m=leaf_plotoptions(m,'layeroffset', 0.1); % distance of secondary layer, generally with cells, to the canvas
m=leaf_plotoptions(m,'gradientoffset', 0.25);%   distance of the polariser gradient arrows from canvas.
m=leaf_plotoptions(m,'highgradcolor',[0,0,0],'lowgradcolor',[0,0,0]); % set colour of polariser gradient arrows for when polariser gradient is high and also when it is low
m=leaf_plotoptions(m,'decorscale',1.2); % size of arrow head
m=leaf_plotoptions(m,'arrowthickness',1.5); % thickness of arrows
m=leaf_plotoptions(m,'sidegrad','AB'); %polariser gradient will be plotted on both sides
%%% END OF USER CODE: INITIALISATION

%%% SECTION 1: ACCESSING MORPHOGENS AND TIME.
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.

    polariser_i = FindMorphogenRole( m, 'POLARISER' );
    P = m.morphogens(:,polariser_i);
    [kapar_i,kapar_p,kapar_a,kapar_l] = getMgenLevels( m, 'KAPAR' );
    [kaper_i,kaper_p,kaper_a,kaper_l] = getMgenLevels( m, 'KAPER' );
    [kbpar_i,kbpar_p,kbpar_a,kbpar_l] = getMgenLevels( m, 'KBPAR' );
    [kbper_i,kbper_p,kbper_a,kbper_l] = getMgenLevels( m, 'KBPER' );
    [knor_i,knor_p,knor_a,knor_l] = getMgenLevels( m, 'KNOR' );
    [strainret_i,strainret_p,strainret_a,strainret_l] = getMgenLevels( m, 'STRAINRET' );
    [arrest_i,arrest_p,arrest_a,arrest_l] = getMgenLevels( m, 'ARREST' );
    [id_base_i,id_base_p,id_base_a,id_base_l] = getMgenLevels( m, 'ID_BASE' );
    [id_tube_i,id_tube_p,id_tube_a,id_tube_l] = getMgenLevels( m, 'ID_TUBE' );
    [id_lobe_i,id_lobe_p,id_lobe_a,id_lobe_l] = getMgenLevels( m, 'ID_LOBE' );
    [id_rim_i,id_rim_p,id_rim_a,id_rim_l] = getMgenLevels( m, 'ID_RIM' );
    [id_med_i,id_med_p,id_med_a,id_med_l] = getMgenLevels( m, 'ID_MED' );
    [id_jun_i,id_jun_p,id_jun_a,id_jun_l] = getMgenLevels( m, 'ID_JUN' );
    [id_sinus_i,id_sinus_p,id_sinus_a,id_sinus_l] = getMgenLevels( m, 'ID_SINUS' );
    [s_jun_i,s_jun_p,s_jun_a,s_jun_l] = getMgenLevels( m, 'S_JUN' );
    [s_med_i,s_med_p,s_med_a,s_med_l] = getMgenLevels( m, 'S_MED' );
    [v_flower_i,v_flower_p,v_flower_a,v_flower_l] = getMgenLevels( m, 'V_FLOWER' );
    [id_tip_i,id_tip_p,id_tip_a,id_tip_l] = getMgenLevels( m, 'ID_TIP' );
    [v_karea_i,v_karea_p,v_karea_a,v_karea_l] = getMgenLevels( m, 'V_KAREA' );
    [v_kaniso_i,v_kaniso_p,v_kaniso_a,v_kaniso_l] = getMgenLevels( m, 'V_KANISO' );
    [id_test_i,id_test_p,id_test_a,id_test_l] = getMgenLevels( m, 'ID_TEST' );
    [id_rep_i,id_rep_p,id_rep_a,id_rep_l] = getMgenLevels( m, 'ID_REP' );
    [id_jcup_i,id_jcup_p,id_jcup_a,id_jcup_l] = getMgenLevels( m, 'ID_JCUP' );
    [id_bcup_i,id_bcup_p,id_bcup_a,id_bcup_l] = getMgenLevels( m, 'ID_BCUP' );
    [id_dorsal_i,id_dorsal_p,id_dorsal_a,id_dorsal_l] = getMgenLevels( m, 'ID_DORSAL' );
    [v_theta_i,v_theta_p,v_theta_a,v_theta_l] = getMgenLevels( m, 'V_THETA' );

% Mesh type: cylinder
%         basecap: 0
%      baseheight: 1
%       baserings: 0
%          centre: 0
%      circumdivs: 80
%          height: 0.02
%      heightdivs: 6
%          layers: 0
%             new: 0
%      randomness: 0
%       thickness: 0
%          topcap: 0
%       topheight: 1
%        toprings: 0
%          xwidth: 0.2
%          ywidth: 0.2

%            Morphogen    Diffusion   Decay   Dilution   Mutant
%            --------------------------------------------------
%                KAPAR         ----    ----       ----     ----
%                KAPER         ----    ----       ----     ----
%                KBPAR         ----    ----       ----     ----
%                KBPER         ----    ----       ----     ----
%                 KNOR         ----    ----       ----     ----
%            POLARISER        0.001    0.01       ----     ----
%            STRAINRET         ----    ----       ----     ----
%               ARREST         ----    ----       ----     ----
%              ID_BASE         ----    ----       ----     ----
%              ID_TUBE         ----    ----       ----     ----
%              ID_LOBE         ----    ----       ----     ----
%               ID_RIM         ----    ----       ----     ----
%               ID_MED         ----    ----       ----     ----
%               ID_JUN         ----    ----       ----     ----
%             ID_SINUS         ----    ----       ----     ----
%                S_JUN         ----    ----       ----     ----
%                S_MED         ----    ----       ----     ----
%             V_FLOWER         ----    ----       ----     ----
%               ID_TIP         ----    ----       ----     ----
%              V_KAREA         ----    ----       ----     ----
%             V_KANISO         ----    ----       ----     ----
%              ID_TEST         ----    ----       ----     ----
%               ID_REP         ----    ----       ----     ----
%              ID_JCUP         ----    ----       ----     ----
%              ID_BCUP         ----    ----       ----     ----
%            ID_DORSAL         ----    ----       ----     ----
%              V_THETA         ----    ----       ----     ----


%%% USER CODE: MORPHOGEN INTERACTIONS



if  realtime == 140
    
    % v_theta is the angle around the z axis.  This is for visualisation
    % purposes and does not model a biological morphogen.
    v_theta_p = atan2( m.nodes(:,2), m.nodes(:,1) )*(0.5/pi);
    
    %specify the base region
    basenodes = m.nodes(:,3)== min(m.nodes(:,3));
    topnodes = m.nodes(:,3)== max(m.nodes(:,3));
    xy = m.nodes(:,[1 2]);
    radiisq = sum(xy.*xy,2);
    maxrsq = max(radiisq);
    
    id_base_p(basenodes) = 1;
    id_base_l = id_base_p;
    m = leaf_fix_vertex( m, 'vertex', basenodes, 'dfs', 'z' );
    
    id_bcup_p = id_base_p;
    
    % create region to be marked with clones
    v_flower_p=ones(size(v_flower_p));
    v_flower_p(basenodes)=0;
    
    % Set id_rim to 1 at every node on the distal edge.
    id_rim_p( topnodes ) = 1 ;
    id_rim_l = id_rim_p .* id_rim_a;
    
    maxlatheight = 0.15;
    
    petalborderv_theta = (1:5)*0.2 - 0.5;
    petalmidv_theta = petalborderv_theta - 0.1;
    v_thetatol = 0.01; %was 0.01
    petalbordernodes = false( size(m.nodes,1), 1 );
    petalmidnodes = false( size(m.nodes,1), 1 );
    for i=1:length(v_theta_p)
        th = v_theta_p(i);
        bth = abs(petalborderv_theta-th);
        if any( (bth < v_thetatol) | (1-bth < v_thetatol) )
            petalbordernodes(i) = true;
        end
        mth = abs(petalmidv_theta-th);
        if any( (mth < v_thetatol) | (1-mth < v_thetatol) )
            petalmidnodes(i) = true;
        end
    end
    id_jun_p(petalbordernodes) = 1;
    id_med_p(:) = zeros(size(id_med_p));
    id_med_p(petalmidnodes) = 1;
    id_jun_l = id_jun_p .* id_jun_a;
    id_med_l = id_med_p .* id_med_a;
    
    lobenodes = m.nodes(:,3) > 0;
    id_lobe_p(lobenodes) = 1;
    % the levels of clamped nodes (vertices) cannot change during
    % simulation
    m.morphogenclamp(lobenodes,id_lobe_i) = 1;
    
    tubenodes = m.nodes(:,3) <= 0;
    id_tube_p(tubenodes) = 1;
    % the levels of clamped nodes (vertices) cannot change during
    % simulation
    m.morphogenclamp(tubenodes,id_tube_i) = 1;
    
    dorsalnodes = m.nodes(:,2) > 0;
    id_dorsal_p(dorsalnodes) = 1;
    % the levels of clamped nodes (vertices) cannot change during
    % simulation
    m.morphogenclamp(dorsalnodes,id_dorsal_i) = 1;
    
    % s_jun is activated by id_jun
    s_jun_p (:)=0;
    s_jun_p (:) = 5* id_jun_l;
    m = leaf_mgen_conductivity( m, 's_jun', 0.00001 );
    m = leaf_mgen_absorption( m, 's_jun', 0.001 );
    
    % s_med is activated by id_jun
    s_med_p (:)=0;
    s_med_p (:) = 5*id_med_l;
    m = leaf_mgen_conductivity( m, 'S_MED', 0.00001 );
    m = leaf_mgen_absorption( m, 'S_MED', 0.001 );
    
    
end


if realtime == 150
    
    m = leaf_mgen_conductivity( m, 's_JUN', 0 );
    m = leaf_mgen_absorption( m, 'S_JUN', 0);
    m = leaf_mgen_conductivity( m, 'S_MED', 0 );
    m = leaf_mgen_absorption( m, 'S_MED', 0 );
    m = leaf_mgen_conductivity( m, 'S_SINUS', 0 );
    m = leaf_mgen_absorption( m, 'S_SINUS', 0 );
    
end



% Code for specific models.
switch modelname
    
    case {'isotropic growth Figure 4B'}
        
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            kapar_p(:) = 0.0165;
            kaper_p(:) = 0.0165;
            kbpar_p(:) = kapar_p;
            kbper_p(:) = kaper_p;
            knor_p(:)  = 0.005;
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
        
    case  'isotropic growth inhibited at junction Figure 4C'   
      
        id_jcup_p= id_jun_p;
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            kapar_p(:) = 0.0165...
                .* inh (100, id_jcup_l);
            
            kaper_p(:) = 0.0165...
                .* inh (100, id_jcup_l);
            
            kbpar_p(:) = kapar_p;
            kbper_p(:) = kaper_p;
            knor_p(:)  = 0.005;
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
        
    case 'isotropic growth inhibited at sinus Figure 4D' % @@model MODEL1
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
        
        id_jcup_p = id_jun_l.* inh(1000, id_tube_l);
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            % Every equation to be formatted should end with an at-at Eqn N comment.
            kapar_p(:) = 0.0165...
                .* inh (100, id_jcup_p);
            kaper_p(:) = 0.0165...
                .* inh (100, id_jcup_p); 
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p;  
            knor_p(:)  = 0.005;  
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
        
    case 'isotropic growth inhibited at sinus and base Figure 4E' % @@model MODEL1
        
        % @@GRN Gene Regulatory Network
     
          id_jcup_p = id_jun_l.* inh(1000, id_tube_l);
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            kapar_p(:) = 0.0165...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);  
            kaper_p(:) = 0.0165...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p); 
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p; 
            knor_p(:)  = 0.005;  
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end

        
    case 'anisotropic growth Figure 4F' % @@model MODEL1
        % @@PRN Polariser Regulatory Network
        
        P(id_base_p ==1)=1; %production of polariser at the base of canvas
        P(id_rim_p==1)=0.1; %degradation of polariser at the edge of the ring of the canvas
        m.morphogenclamp((id_rim_p==1)|(id_base_p==1),polariser_i) = 1;%stable polariser gradient by clamping the production and degradation of the polariser
        m = leaf_mgen_conductivity( m, 'Polariser', 0.001);% diffusion constant for polariser
        m = leaf_mgen_dilution( m, 'Polariser', false );% it will not dilute with growth
        m = leaf_mgen_absorption( m, 'Polariser', 0.01);  % it will not decay everywhere
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
        
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350.
            kapar_p(:) = 0.022;
            kaper_p(:) = 0.011;
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p;  
            knor_p(:)  = 0.005;  
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
        
    case 'anisotropic growth inhibited at junction Figure 4G'   % @@model MODEL1
        % @@PRN Polariser Regulatory Network
        
        P(id_base_p ==1)=1; %production of polariser at the base of canvas
        P(id_rim_p==1)=0.1; %degradation of polariser at the edge of the ring of the canvas
        m.morphogenclamp((id_rim_p==1)|(id_base_p==1),polariser_i) = 1;%stable polariser gradient by clamping the production and degradation of the polariser
        m = leaf_mgen_conductivity( m, 'Polariser', 0.001);% diffusion constant for polariser
        m = leaf_mgen_dilution( m, 'Polariser', false );% it will not dilute with growth
        m = leaf_mgen_absorption( m, 'Polariser', 0.01);  % it will not decay everywhere
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
       
        
         id_jcup_p= id_jun_p;
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            % Every equation to be formatted should end with an at-at Eqn N comment.
            kapar_p(:) = 0.022...
                .* inh (100, id_jcup_p);
            kaper_p(:) = 0.011...
                .* inh (100, id_jcup_p);    
            kbpar_p(:) = kapar_p;  % @@ Eqn xx
            kbper_p(:) = kaper_p;  % @@ Eqn xx
            knor_p(:)  = 0.01;  % @@ Eqn xx increase from 0.005
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
    case 'anisotropic growth inhibited at sinus Figure 4H'   % @@model MODEL1
        % @@PRN Polariser Regulatory Network
        
        P(id_base_p ==1)=1; %production of polariser at the base of canvas
        P(id_rim_p==1)=0.1; %degradation of polariser at the edge of the ring of the canvas
        m.morphogenclamp((id_rim_p==1)|(id_base_p==1),polariser_i) = 1;%stable polariser gradient by clamping the production and degradation of the polariser
        m = leaf_mgen_conductivity( m, 'Polariser', 0.001);% diffusion constant for polariser
        m = leaf_mgen_dilution( m, 'Polariser', false );% it will not dilute with growth
        m = leaf_mgen_absorption( m, 'Polariser', 0.01);  % it will not decay everywhere
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
      
          id_jcup_p = id_jun_l.* inh(1000, id_tube_l);
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            % Every equation to be formatted should end with an at-at Eqn N comment.
            kapar_p(:) = 0.022...
                .* inh (100, id_jcup_p);
            kaper_p(:) = 0.011...
                .* inh (100, id_jcup_p);
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p;  
            knor_p(:)  = 0.005;  
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
    case 'anisotropic growth inhibited at sinus and base Figure 4I' % @@model MODEL1
        % @@PRN Polariser Regulatory Network
        
        P(id_base_p ==1)=1; %production of polariser at the base of canvas
        P(id_rim_p==1)=0.1; %degradation of polariser at the edge of the ring of the canvas
        m.morphogenclamp((id_rim_p==1)|(id_base_p==1),polariser_i) = 1;%stable polariser gradient by clamping the production and degradation of the polariser
        m = leaf_mgen_conductivity( m, 'Polariser', 0.001);% diffusion constant for polariser
        m = leaf_mgen_dilution( m, 'Polariser', false );% it will not dilute with growth
        m = leaf_mgen_absorption( m, 'Polariser', 0.01);  % it will not decay everywhere
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
        
        
         id_jcup_p = id_jun_l.* inh(1000, id_tube_l);
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            % Every equation to be formatted should end with an at-at Eqn N comment.
            kapar_p(:) = 0.022...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);
            kaper_p(:) = 0.011...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);  
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p;  
            knor_p(:)  = 0.005; 
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
         case 'growth promoted by med Figure 4K' % @@model MODEL1
        % @@PRN Polariser Regulatory Network
        
        P(id_base_p ==1)=1; %production of polariser at the base of canvas
        P(id_rim_p==1)=0.1; %degradation of polariser at the edge of the ring of the canvas
        m.morphogenclamp((id_rim_p==1)|(id_base_p==1),polariser_i) = 1;%stable polariser gradient by clamping the production and degradation of the polariser
        m = leaf_mgen_conductivity( m, 'Polariser', 0.001);% diffusion constant for polariser
        m = leaf_mgen_dilution( m, 'Polariser', false );% it will not dilute with growth
        m = leaf_mgen_absorption( m, 'Polariser', 0.01);  % it will not decay everywhere
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
       
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            % Every equation to be formatted should end with an at-at Eqn N comment.
            kapar_p(:) = 0.022...was 0.011
                .*pro(0.3, s_med_l)...
                .* inh (0.5, id_base_p);
               
            kaper_p(:) = 0.011...
                .*pro(0, s_med_l)...
                .* inh (0.5, id_base_p);
                
                  
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p;  
            knor_p(:)  = 0.005; 
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
        case 'growth promoted by med and inhibited at junction Figure 4L' % @@model MODEL1
        % @@PRN Polariser Regulatory Network
        
        P(id_base_p ==1)=1; %production of polariser at the base of canvas
        P(id_rim_p==1)=0.1; %degradation of polariser at the edge of the ring of the canvas
        m.morphogenclamp((id_rim_p==1)|(id_base_p==1),polariser_i) = 1;%stable polariser gradient by clamping the production and degradation of the polariser
        m = leaf_mgen_conductivity( m, 'Polariser', 0.001);% diffusion constant for polariser
        m = leaf_mgen_dilution( m, 'Polariser', false );% it will not dilute with growth
        m = leaf_mgen_absorption( m, 'Polariser', 0.01);  % it will not decay everywhere
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
        
         id_jcup_p= id_jun_p;
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            % Every equation to be formatted should end with an at-at Eqn N comment.
            kapar_p(:) = 0.022...was 0.011
                .*pro(0.3, s_med_l)...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);
       
            kaper_p(:) = 0.011...
                .*pro(0, s_med_l)...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);
                
                  
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p;  
            knor_p(:)  = 0.005; 
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
        case 'growth promoted by med and inhibited at sinus Figure 4M' % @@model MODEL1
        % @@PRN Polariser Regulatory Network
        
        P(id_base_p ==1)=1; %production of polariser at the base of canvas
        P(id_rim_p==1)=0.1; %degradation of polariser at the edge of the ring of the canvas
        m.morphogenclamp((id_rim_p==1)|(id_base_p==1),polariser_i) = 1;%stable polariser gradient by clamping the production and degradation of the polariser
        m = leaf_mgen_conductivity( m, 'Polariser', 0.001);% diffusion constant for polariser
        m = leaf_mgen_dilution( m, 'Polariser', false );% it will not dilute with growth
        m = leaf_mgen_absorption( m, 'Polariser', 0.01);  % it will not decay everywhere
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
        
        id_jcup_p = id_jun_l.* inh(1000, id_tube_l);
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            % Every equation to be formatted should end with an at-at Eqn N comment.
            kapar_p(:) = 0.022...
                .*pro(0.3, s_med_l)...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);
               
            kaper_p(:) = 0.011...
                .*pro(0, s_med_l)...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);
                
                  
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p;  
            knor_p(:)  = 0.005; 
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
        

        
        case  'growth promoted by med and selectively inhibited at sinus Figure 4N' % @@model MODEL1
        % @@PRN Polariser Regulatory Network
        
        P(id_base_p ==1)=1; %production of polariser at the base of canvas
        P(id_rim_p==1)=0.1; %degradation of polariser at the edge of the ring of the canvas
        m.morphogenclamp((id_rim_p==1)|(id_base_p==1),polariser_i) = 1;%stable polariser gradient by clamping the production and degradation of the polariser
        m = leaf_mgen_conductivity( m, 'Polariser', 0.001);% diffusion constant for polariser
        m = leaf_mgen_dilution( m, 'Polariser', false );% it will not dilute with growth
        m = leaf_mgen_absorption( m, 'Polariser', 0.01);  % it will not decay everywhere
        
        % @@GRN Gene Regulatory Network
        % Every equation to be formatted should end with an at-at Eqn N comment.
        
        id_jcup_p = id_jun_l.* inh(1000, id_tube_l).*id_dorsal_p;
        
        % @@KRN Growth Regulatory Network
        if realtime == 160 && realtime<350
            % Every equation to be formatted should end with an at-at Eqn N comment.
            kapar_p(:) = 0.022...
                .*pro(0.3, s_med_l)...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);
               
            kaper_p(:) = 0.011...
                .*pro(0, s_med_l)...
                .* inh (100, id_jcup_p)...
                .* inh (0.5, id_base_p);
                
                  
            kbpar_p(:) = kapar_p;  
            kbper_p(:) = kaper_p;  
            knor_p(:)  = 0.005; 
            v_karea_p = kapar_p +kaper_p;
            v_kaniso_p =  log (kbpar_p./kbper_p);
        end
 otherwise
        % If this happens, maybe you forgot a model.

end
%%% END OF USER CODE: MORPHOGEN INTERACTIONS

%%% SECTION 3: INSTALLING MODIFIED VALUES BACK INTO MESH STRUCTURE
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.
    m.morphogens(:,polariser_i) = P;
    m.morphogens(:,kapar_i) = kapar_p;
    m.morphogens(:,kaper_i) = kaper_p;
    m.morphogens(:,kbpar_i) = kbpar_p;
    m.morphogens(:,kbper_i) = kbper_p;
    m.morphogens(:,knor_i) = knor_p;
    m.morphogens(:,strainret_i) = strainret_p;
    m.morphogens(:,arrest_i) = arrest_p;
    m.morphogens(:,id_base_i) = id_base_p;
    m.morphogens(:,id_tube_i) = id_tube_p;
    m.morphogens(:,id_lobe_i) = id_lobe_p;
    m.morphogens(:,id_rim_i) = id_rim_p;
    m.morphogens(:,id_med_i) = id_med_p;
    m.morphogens(:,id_jun_i) = id_jun_p;
    m.morphogens(:,id_sinus_i) = id_sinus_p;
    m.morphogens(:,s_jun_i) = s_jun_p;
    m.morphogens(:,s_med_i) = s_med_p;
    m.morphogens(:,v_flower_i) = v_flower_p;
    m.morphogens(:,id_tip_i) = id_tip_p;
    m.morphogens(:,v_karea_i) = v_karea_p;
    m.morphogens(:,v_kaniso_i) = v_kaniso_p;
    m.morphogens(:,id_test_i) = id_test_p;
    m.morphogens(:,id_rep_i) = id_rep_p;
    m.morphogens(:,id_jcup_i) = id_jcup_p;
    m.morphogens(:,id_bcup_i) = id_bcup_p;
    m.morphogens(:,id_dorsal_i) = id_dorsal_p;
    m.morphogens(:,v_theta_i) = v_theta_p;

%%% USER CODE: FINALISATION

% In this section you may modify the mesh in any way whatsoever.

% If needed force FE to subdivide (increase number FE's) here
% if realtime==280+dt
% m = leaf_subdivide( m, 'morphogen','id_vent',...
%       'min',0.5,'max',1,...
%       'mode','mid','levels','all');
% end
% Cut the mesh along the seams (see above)
% if m.userdata.CutOpen==1
%    m=leaf_dissect(m);
%    m.userdata.CutOpen=2;
%    Relax accumulated stresses slowly i.e. 0.95 to 0.999
%    m = leaf_setproperty( m, 'freezing', 0.999 );
% end
%%% END OF USER CODE: FINALISATION

end


%%% USER CODE: SUBFUNCTIONS

function m = local_setproperties( m )
% This function is called at time zero in the INITIALISATION section of the
% interaction function.  It provides commands to set each of the properties
% that are contained in m.globalProps.  Uncomment whichever ones you would
% like to set yourself, and put in whatever value you want.
%
% Some of these properties are for internal use only and should never be
% set by the user.  At some point these will be moved into a different
% component of m, but for the present, just don't change anything unless
% you know what it is you're changing.

%    m = leaf_setproperty( m, 'trinodesvalid', true );
%    m = leaf_setproperty( m, 'prismnodesvalid', true );
%    m = leaf_setproperty( m, 'thresholdsq', 0.003053 );
%    m = leaf_setproperty( m, 'lengthscale', 0.199662 );
%    m = leaf_setproperty( m, 'initialArea', 0.068817 );
%    m = leaf_setproperty( m, 'bendunitlength', 0.262330 );
%    m = leaf_setproperty( m, 'thicknessRelative', 0.100000 );
%    m = leaf_setproperty( m, 'thicknessArea', 1.000000 );
%    m = leaf_setproperty( m, 'hybridMesh', false );
%    m = leaf_setproperty( m, 'thicknessMode', 'physical' );
%    m = leaf_setproperty( m, 'activeGrowth', 1.000000 );
%    m = leaf_setproperty( m, 'displayedGrowth', 1.000000 );
%    m = leaf_setproperty( m, 'displayedMulti', [] );
%    m = leaf_setproperty( m, 'allowNegativeGrowth', true );
%    m = leaf_setproperty( m, 'usePrevDispAsEstimate', true );
%    m = leaf_setproperty( m, 'perturbInitGrowthEstimate', 0.000010 );
%    m = leaf_setproperty( m, 'perturbRelGrowthEstimate', 0.010000 );
%    m = leaf_setproperty( m, 'perturbDiffusionEstimate', 0.000100 );
%    m = leaf_setproperty( m, 'resetRand', false );
%    m = leaf_setproperty( m, 'mingradient', 0.000000 );
%    m = leaf_setproperty( m, 'relativepolgrad', false );
%    m = leaf_setproperty( m, 'usefrozengradient', true );
%    m = leaf_setproperty( m, 'userpolarisation', false );
%    m = leaf_setproperty( m, 'twosidedpolarisation', false );
%    m = leaf_setproperty( m, 'splitmargin', 1.400000 );
%    m = leaf_setproperty( m, 'splitmorphogen', '' );
%    m = leaf_setproperty( m, 'thresholdmgen', 0.500000 );
%    m = leaf_setproperty( m, 'bulkmodulus', 1.000000 );
%    m = leaf_setproperty( m, 'unitbulkmodulus', true );
%    m = leaf_setproperty( m, 'poissonsRatio', 0.300000 );
%    m = leaf_setproperty( m, 'starttime', 0.000000 );
%    m = leaf_setproperty( m, 'timestep', 0.010000 );
%    m = leaf_setproperty( m, 'timeunitname', '' );
%    m = leaf_setproperty( m, 'distunitname', 'mm' );
%    m = leaf_setproperty( m, 'scalebarvalue', 0.000000 );
%    m = leaf_setproperty( m, 'validateMesh', true );
%    m = leaf_setproperty( m, 'rectifyverticals', false );
%    m = leaf_setproperty( m, 'allowSplitLongFEM', true );
%    m = leaf_setproperty( m, 'allowSplitThinFEM', false );
%    m = leaf_setproperty( m, 'splitthinness', 10.000000 );
%    m = leaf_setproperty( m, 'longSplitThresholdPower', 0.000000 );
%    m = leaf_setproperty( m, 'allowSplitBentFEM', false );
%    m = leaf_setproperty( m, 'allowSplitBio', true );
%    m = leaf_setproperty( m, 'allowFlipEdges', false );
%    m = leaf_setproperty( m, 'allowElideEdges', true );
%    m = leaf_setproperty( m, 'mincellangle', 0.200000 );
%    m = leaf_setproperty( m, 'mincellrelarea', 0.040000 );
%    m = leaf_setproperty( m, 'alwaysFlat', 0.000000 );
%    m = leaf_setproperty( m, 'flattenforceconvex', true );
%    m = leaf_setproperty( m, 'flatten', false );
%    m = leaf_setproperty( m, 'flattenratio', 1.000000 );
%    m = leaf_setproperty( m, 'useGrowthTensors', false );
%    m = leaf_setproperty( m, 'useMorphogens', true );
%    m = leaf_setproperty( m, 'plasticGrowth', false );
%    m = leaf_setproperty( m, 'maxFEcells', 0 );
%    m = leaf_setproperty( m, 'inittotalcells', 0 );
%    m = leaf_setproperty( m, 'bioApresplitproc', '' );
%    m = leaf_setproperty( m, 'bioApostsplitproc', '' );
%    m = leaf_setproperty( m, 'maxBioAcells', 0 );
%    m = leaf_setproperty( m, 'biosplitarea', 0.000000 );
%    m = leaf_setproperty( m, 'biosplitarrestmgen', 'ARREST' );
%    m = leaf_setproperty( m, 'biosplitarrestmgenthreshold', 0.990000 );
%    m = leaf_setproperty( m, 'bioMinEdgeLength', 0.000000 );
%    m = leaf_setproperty( m, 'bioSpacePullInRatio', 0.100000 );
%    m = leaf_setproperty( m, 'colors', (6 values) );
%    m = leaf_setproperty( m, 'colorvariation', 0.050000 );
%    m = leaf_setproperty( m, 'colorparams', (12 values) );
%    m = leaf_setproperty( m, 'biocolormode', 'auto' );
%    m = leaf_setproperty( m, 'userpostiterateproc', [] );
%    m = leaf_setproperty( m, 'canceldrift', false );
%    m = leaf_setproperty( m, 'mgen_interaction', '' );
%    m = leaf_setproperty( m, 'mgen_interactionName', 'gpt_petal_primordia_ring_20160616' );
%    m = leaf_setproperty( m, 'allowInteraction', true );
%    m = leaf_setproperty( m, 'interactionValid', true );
%    m = leaf_setproperty( m, 'gaussInfo', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'D', (36 values) );
%    m = leaf_setproperty( m, 'C', (36 values) );
%    m = leaf_setproperty( m, 'G', (6 values) );
%    m = leaf_setproperty( m, 'solver', 'cgs' );
%    m = leaf_setproperty( m, 'solverprecision', 'double' );
%    m = leaf_setproperty( m, 'solvertolerance', 0.001000 );
%    m = leaf_setproperty( m, 'solvertolerancemethod', 'max' );
%    m = leaf_setproperty( m, 'diffusiontolerance', 0.000010 );
%    m = leaf_setproperty( m, 'allowsparse', true );
%    m = leaf_setproperty( m, 'maxsolvetime', 1000.000000 );
%    m = leaf_setproperty( m, 'cgiters', 0 );
%    m = leaf_setproperty( m, 'simsteps', 0 );
%    m = leaf_setproperty( m, 'stepsperrender', 0 );
%    m = leaf_setproperty( m, 'growthEnabled', true );
%    m = leaf_setproperty( m, 'diffusionEnabled', true );
%    m = leaf_setproperty( m, 'flashmovie', false );
%    m = leaf_setproperty( m, 'makemovie', false );
%    m = leaf_setproperty( m, 'moviefile', '' );
%    m = leaf_setproperty( m, 'codec', 'Motion JPEG AVI' );
%    m = leaf_setproperty( m, 'autonamemovie', true );
%    m = leaf_setproperty( m, 'overwritemovie', false );
%    m = leaf_setproperty( m, 'framesize', [] );
%    m = leaf_setproperty( m, 'mov', [] );
%    m = leaf_setproperty( m, 'boingNeeded', false );
%    m = leaf_setproperty( m, 'defaultinterp', 'min' );
%    m = leaf_setproperty( m, 'readonly', false );
%    m = leaf_setproperty( m, 'projectdir', 'C:\Users\xana\Desktop\Modelling' );
%    m = leaf_setproperty( m, 'modelname', 'GPT_petal_primordia_ring_20160616' );
%    m = leaf_setproperty( m, 'allowsave', true );
%    m = leaf_setproperty( m, 'addedToPath', false );
%    m = leaf_setproperty( m, 'bendsplit', 0.300000 );
%    m = leaf_setproperty( m, 'usepolfreezebc', false );
%    m = leaf_setproperty( m, 'dorsaltop', true );
%    m = leaf_setproperty( m, 'defaultazimuth', -45.000000 );
%    m = leaf_setproperty( m, 'defaultelevation', 33.750000 );
%    m = leaf_setproperty( m, 'defaultroll', 0.000000 );
%    m = leaf_setproperty( m, 'defaultViewParams', (unknown type ''struct'') );
%    m = leaf_setproperty( m, 'comment', '' );
%    m = leaf_setproperty( m, 'legendTemplate', '%T: %q\n%m' );
%    m = leaf_setproperty( m, 'bioAsplitcells', true );
%    m = leaf_setproperty( m, 'bioApullin', 0.142857 );
%    m = leaf_setproperty( m, 'bioAfakepull', 0.202073 );
%    m = leaf_setproperty( m, 'viewrotationstart', -45.000000 );
%    m = leaf_setproperty( m, 'viewrotationperiod', 0.000000 );
%    m = leaf_setproperty( m, 'interactive', false );
%    m = leaf_setproperty( m, 'coderevision', 0 );
%    m = leaf_setproperty( m, 'coderevisiondate', '' );
%    m = leaf_setproperty( m, 'modelrevision', 0 );
%    m = leaf_setproperty( m, 'modelrevisiondate', '' );
%    m = leaf_setproperty( m, 'savedrunname', '' );
%    m = leaf_setproperty( m, 'savedrundesc', '' );
%    m = leaf_setproperty( m, 'vxgrad', (108 values) );
end

% Here you may write any functions of your own, that you want to call from
% the interaction function, but never need to call from outside it.
% Remember that they do not have access to any variables except those
% that you pass as parameters, and cannot change anything except by
% returning new values as results.
% Whichever section they are called from, they must respect the same
% restrictions on what modifications they are allowed to make to the mesh.

% For example:

% function m = do_something( m )
%   % Change m in some way.
% end

% Call it from the main body of the interaction function like this:
%       m = do_something( m );
