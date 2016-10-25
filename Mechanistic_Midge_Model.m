function Mechanistic_Midge_Model
%MECHANISTIC_MIDGE_MODEL Midge model using example data and plotted
%   Note that this model uses example data, but daily field data may be
%   used.  Here we use a leap year number of days to coincide with our
%   field data.
%
%   It is assumed that all data files are in the same folder as this
%   script.
%
%   This was programmed by Dr Steven White, CEH Wallingford (25/10/16).

%% Load datasets
tempvar=load('photoperiod.mat'); %photoperiod data
photoperiod=tempvar.photoperiod;

tempvar=load('soiltemps.mat'); %load the soil temps
soiltemps=tempvar.soiltemps;

tempvar=load('airtemps.mat'); %load the temperature data.
airtemps=tempvar.airtemps;

ddp=[1.33292201187087e-08; 2.67591445332698; 3.18698597432782e-07; ...
    5.69639631283444; 8.64223690531383e-05; 46986.5044460894; ...
    49065.5221025845]; %Define site specific model parameters

%% Run the model and plot

[totalmodeloutput,ages]=model(ddp); %call the main model function below

figure(1)
plot(totalmodeloutput,'r','LineWidth',2)
xlim([1 366])
xlabel('Days since 1st Jan 2008','Fontsize',12)
ylabel('Adult Abundance','Fontsize',12)
title('Example Model Output','FontSize',12)

figure(2)
hold on
for index=1:366
    agetp=ages{index,1};
    agetp2=unique(agetp);
    plot(index*ones(size(agetp2)),agetp2,'k.')
end
axis([1 366 0 1.1])
xlabel('Days since 1st Jan 2008','Fontsize',12)
ylabel('Pre-Adult Age','Fontsize',12)
title('Example Model Output','Fontsize',12)



    function [modeloutput,ages]=model(ddparams1)
        % Define the initial cohorts of midges. Here we assume that
        % initially all midges are pre-pupal.
        
        % Define the initial number of cohorts.
        ncoh=50;
        % Define the number of eggs in each cohort (comes from data). There
        % are no adults
        eggnum=20*normrnd(49.7,11.1,ncoh,1); %number of adults*clutch size
        nullnum=[];
        parnum=[];
        % Next we define the pre-adult age of each cohort.  To do this, we
        % define age on the interval [0 1], where 0 indicates new eggs,
        % 0.91 indicates new pupae and 1 indicates adult emergence from the
        % pupal stage.  The ages come from data.
        eggage=zeros(ncoh,1); %Initialise the eggages.  There are no adults.
        nullage=[];
        parage=[];
        for i=1:ncoh
            agetemporary=normrnd(0.76,0.2); %pick an age
            while agetemporary>=0.91 %check to see if midges are pre-pupal
                agetemporary=normrnd(0.76,0.1); %if not pre-pupal try a new age
            end
            eggage(i)=agetemporary; %set age
        end
        clear agetemporary
        
        
        % Number of days in year.
        timeend=length(airtemps);
        
        % Do the main model run
        
        % Run the model once to recalculate the initial pre-adult numbers
        % and ages
        [~,~,~,ICsegg,~,~,~,~,~,~]=yearrun(airtemps,timeend,eggnum,nullnum,parnum,eggage,nullage,parage);
        % do a second run for luck ;-)
        [~,~,~,ICsegg,~,~,~,~,~,~]=yearrun(airtemps,timeend,ICsegg(:,1),[],[],ICsegg(:,2),[],[]);
        % Do the proper run
        [~,nullvec,parvec,~,~,~,~,~,~,ages]=yearrun(airtemps,timeend,ICsegg(:,1),[],[],ICsegg(:,2),[],[]);
        modeloutput=ddparams1(3)*(nullvec+parvec); %adjust output by efficiency coefficient
        
        %% Yearly model run - this is the main function
        function [eggvec,nullvec,parvec,ICsegg,ICsnull,ICspar,newegg,newnull,newpar,ages]=yearrun(temp,timeend,eggnum,nullnum,parnum,eggage,nullage,parage)
            
            % Define a load of temporary storage vectors (used for moving
            % cohorts between clases).
            nullfromeggage=[];
            nullfromeggnum=[];
            parfromnullage=[];
            parfromnullnum=[];
            eggfromnullage=[];
            eggfromnullnum=[];
            parfromparage=[];
            parfromparnum=[];
            eggfromparage=[];
            eggfromparnum=[];
            
            % Define some storage vectors
            eggvec=zeros(timeend,1);
            eggvec(1)=sum(eggnum);
            nullvec=zeros(timeend,1);
            nullvec(1)=sum(nullnum);
            parvec=zeros(timeend,1);
            parvec(1)=sum(parnum);
            
            newegg=[];
            newnull=[];
            newpar=[];
            
            track=zeros(timeend+1,2);
            track(1,1)=max(eggage);
            
            ages=cell(timeend,3);
            
            marker=1;
            
            for t=1:timeend
                
                if marker==1 && track(t,1)<1
                    [track(t+1,1),~]=eggdev(track(t,1),0,t,ddparams1);
                elseif marker==1 && track(t,1)>1
                    track(t+1,1)=0;
                    marker=2;
                elseif marker==2
                    marker=3;
                elseif marker==3 && track(t,2)<1
                    [track(t+1,2),~]=adultdev(track(t,2),0,t);
                elseif marker==3 && track(t,2)>1
                    track(t+1,2)=0;
                    marker=1;
                end
                
                %Check to see if any of populations have gone extinct
                if nnz(parnum<=0.001)>0 %Check to see if the parous adults have gone extinct
                    parage=parage(parnum>0.001); %Remove the ages of parous adults that have gone extinct
                    parnum=parnum(parnum>0.001); %Remove the populations parous adults that have gone extinct
                end
                if nnz(nullnum<=0.001)>0 %Check to see if any of the nulliparous adults have gone extinct
                    nullage=nullage(nullnum>0.001); %Remove the ages of nulliparous adults that have gone extinct
                    nullnum=nullnum(nullnum>0.001); %Remove the populations nulliparous adults that have gone extinct
                end
                
                % Evolve the cohorts
                
                % Eggs (pre-adults)
                if nnz(eggage>=1)==0  && ~isempty(eggage) %Some eggs are available to develop
                    [eggage,eggnum]=eggdev(eggage,eggnum,t,ddparams1); %Develop the pre-adults
                    nullfromeggage=[]; %No nulliparous adults come from eggs
                    nullfromeggnum=[]; %No nulliparous adults come from eggs
                elseif nnz(eggage>=1)>0  && ~isempty(eggage) %Some eggs have developed
                    nullfromeggage=zeros(size(eggage(eggage>=1))); %Make new nulliparous adult ages
                    nullfromeggnum=eggnum(eggage>=1); %Add the pre-adult numbers to the nulliparous adults
                    eggnum=eggnum(eggage<1); %Remove the aged pre-adults that have developed from the pre-adult numbers
                    eggage=eggage(eggage<1); %Remove the aged pre-adults that have developed from the pre-adult ages
                    if ~isempty(eggage)
                        [eggage,eggnum]=eggdev(eggage,eggnum,t,ddparams1); %Develop the pre-adults
                    end
                else
                    %There are no pre-adults
                    nullfromeggage=[]; %No nulliparous adults come from eggs
                    nullfromeggnum=[]; %No nulliparous adults come from eggs
                end
                
                % nulliparous Adults
                if nnz(nullage>=1)==0  && ~isempty(nullage) %Some nulliparous adults are available to develop but don't lay eggs
                    [nullage,nullnum]=adultdev(nullage,nullnum,t); %Develop the nulliparous adults
                    parfromnullage=[]; %No parous adults come from nulliparous adults
                    parfromnullnum=[]; %No parous adults come from nulliparous adults
                    eggfromnullage=[]; %No eggs come from nulliparous adults
                    eggfromnullnum=[]; %No eggs come from nulliparous adults
                elseif nnz(nullage>=1)>0  && ~isempty(nullage) %Some nulliparous have developed and lay eggs
                    parfromnullage=zeros(size(nullage(nullage>=1))); %Make new parous adult ages
                    parfromnullnum=nullnum(nullage>=1); %Add the nulliparous adult numbers to the parous adults
                    eggfromnullage=zeros(nnz(nullage>=1),1); %Add ages of newly laid eggs
                    if temp(t)>12
                        eggfromnullnum=nullnum(nullage>=1).*normrnd(49.7,11.1,nnz(nullage>=1),1);%49.7.*ones(nnz(nullage>=1),1); %Add newly laid eggs if warm
                    else
                        eggfromnullnum=nullnum(nullage>=1).*5.*ones(nnz(nullage>=1),1); %Add newly laid eggs if cold
                    end
                    nullnum=nullnum(nullage<1); %Remove the aged nulliparous adults that have developed from the nulliparous adult numbers
                    nullage=nullage(nullage<1); %Remove the aged nulliparous adults that have developed from the nulliparous adult ages
                    if ~isempty(nullage)
                        [nullage,nullnum]=adultdev(nullage,nullnum,t); %Develop the nulliparous adults
                    end
                else
                    %There are no nulliparours adults
                    parfromnullage=[]; %No parous adults come from nulliparous adults
                    parfromnullnum=[]; %No parous adults come from nulliparous adults
                    eggfromnullage=[]; %No eggs come from nulliparous adults
                    eggfromnullnum=[]; %No eggs come from nulliparous adults
                end
                
                % Parous Adults
                if nnz(parage>=1)==0  && ~isempty(parage) %Some parous adults are available to develop but don't lay eggs
                    [parage,parnum]=adultdev(parage,parnum,t); %Develop the parous adults
                    parfromparage=[]; %No parous adults come from parous adults
                    parfromparnum=[]; %No parous adults come from parous adults
                    eggfromparage=[]; %No eggs come from parous adults
                    eggfromparnum=[]; %No eggs come from parous adults
                elseif nnz(parage>=1)>0  && ~isempty(parage) %Some parous have developed and lay eggs
                    parfromparage=zeros(size(parage(parage>=1))); %Make new parous adult ages
                    parfromparnum=parnum(parage>=1); %Add the parous adult numbers to the parous adults
                    eggfromparage=zeros(nnz(parage>=1),1); %Add ages of newly laid eggs
                    if temp(t)>12
                        eggfromparnum=parnum(parage>=1).*normrnd(49.7,11.1,nnz(parage>=1),1);%49.7.*ones(nnz(parage>=1),1); %Add newly laid eggs if warm
                    else
                        eggfromparnum=parnum(parage>=1).*5.*ones(nnz(parage>=1),1); %Add newly laid eggs if cold
                    end
                    parnum=parnum(parage<1); %Remove the aged parous adults that have developed from the parous adult numbers
                    parage=parage(parage<1); %Remove the aged parous adults that have developed from the parous adult ages
                    if ~isempty(parage)
                        [parage,parnum]=adultdev(parage,parnum,t); %Develop the parous adults
                    end
                else
                    %There are no parous adults
                    parfromparage=[]; %No parous adults come from parous adults
                    parfromparnum=[]; %No parous adults come from parous adults
                    eggfromparage=[]; %No eggs come from parous adults
                    eggfromparnum=[]; %No eggs come from parous adults
                end
                
                %Add in all the new eggs, nulliparous adults and parous
                %adults that have been created in this time-step, if any.
                %If none, then only a blank is added!  Do this here to
                %avoid new individuals evolving during a time-step.
                
                eggage=[eggage; eggfromnullage; eggfromparage]; %Egg ages
                eggnum=[eggnum; eggfromnullnum; eggfromparnum]; %Egg numbers
                nullage=[nullage; nullfromeggage]; %nulliparous adult ages
                nullnum=[nullnum; nullfromeggnum]; %nulliparous adult ages
                parage=[parage; parfromnullage; parfromparage]; %Parous adult ages
                parnum=[parnum; parfromnullnum; parfromparnum]; %Parous adult ages
                
                newegg=[newegg; sum(eggfromnullnum)+sum(eggfromparnum)]; %Store new pre-adults
                newnull=[newnull; sum(nullfromeggnum)]; %Store new nulliparous adults
                newpar=[newpar; sum(parfromnullnum)+sum(parfromparnum)]; %Sore new parous adults
                
                % Store the total abundances for the day - note the IC is
                % not kept
                eggvec(t)=sum(eggnum);
                nullvec(t)=sum(nullnum);
                parvec(t)=sum(parnum);
                
                % Store the ages
                ages{t,1}=eggage;
                ages{t,2}=nullage;
                ages{t,3}=parage;
            end
            %Store the last numbers and ages (used for running the function
            %multiple times)
            ICsegg=[eggnum, eggage];
            ICsnull=[nullnum, nullage];
            ICspar=[parnum, parage];
            
        end
        
        
        
        %% Development functions defined in the main function
        
        % Egg development and survival
        function [eggage,eggnum]=eggdev(eggage,eggnum,t,ddparams1)
            %Development
            if soiltemps(t)>=0 %Is it warm enough to develop?
                %Check to see if not in diapause period
                if mod(t,366)>find(photoperiod>ddparams1(end-1),1,'first') && ...
                        mod(t,366)<find(photoperiod>ddparams1(end),1,'last')
                    eggage=eggage+exp(0.1011*soiltemps(t)-5.9177); %If yes then age
                else
                    %only young pre-adults develop
                    eggage(eggage<0.91)=eggage(eggage<0.91)+exp(0.1011*soiltemps(t)-5.9177);
                end
                %Otherwise don't age
            end
            %Mortalities
            eggnum=eggnum*eggsurv(t)*DD(eggnum,t,ddparams1);
            
            % Egg survival
            function alpha=eggsurv(t)
                b1=0.0009;
                c1=0.9629;
                b2=-0.0039;
                c2=1.0802;
                if soiltemps(t)<5
                    alpha=1;
                elseif soiltemps(t)>=5 && soiltemps(t)<12.5
                    alpha=0.974;
                elseif soiltemps(t)>=12.5 && soiltemps(t)<24.5
                    alpha=b1*soiltemps(t)+c1;
                elseif soiltemps(t)>=24.5 && soiltemps(t)<35
                    alpha=b2*soiltemps(t)+c2;
                else
                    alpha=0.944;
                end
            end
            
            % Larval density-depend - note egg/pupal stages incur DD
            % Maynard Smith and Slatkin 1973 Nt+1=r*Nt/(1+(a*Nt)^b)
            function dd=DD(eggnum,t,ddparams1)
                a=ddparams1(1);
                b=ddparams1(2);
                if soiltemps(t)<5
                    dd=1;
                else
                    dd=1/(1+(a*sum(eggnum))^b);
                end
            end
            
        end
        
        % Adult development and survival
        function [adultage,adultnum]=adultdev(adultage,adultnum,t)
            %Development
            if airtemps(t)>=1 %Is it warm enough to develop?
                func=-1.98+0.07217*airtemps(t)+2516.65/(airtemps(t)^2);
                adultage=adultage+1/func; %If yes then age
                %Otherwise don't age
            end
            %Mortalities
            adultnum=adultnum*adultsurv(t);
            
            % Adult survival
            function alpha=adultsurv(t)
                u=0.00001;
                v=-0.001;
                w=0.0187;
                z=0.8924;
                if airtemps(t)<=2
                    alpha=0.2;
                elseif airtemps(t)>2 && airtemps(t)<10
                    alpha=0.9894;
                else
                    alpha=u*airtemps(t)^3+v*airtemps(t)^2+w*airtemps(t)+z;
                end
            end
            
        end
        %% End of Function
    end
end

