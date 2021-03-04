function [returnVar,msg] = RemovePC1_datFile(basepath)

%Based on ReRefFilFile, removes first PC of DAT file
% Eliezyer de Oliveira 2020

cd(basepath)

datFiles = dir('*.dat');
%all dat files you don't want to run this script on, add more as you need.
dats2notuse = {'analogin.dat','auxiliary.dat','digitalin.dat','supply.dat','time.dat'};
%we only use the dat file recording the amplifier
dat2use = ~ismember({datFiles.name},dats2notuse);
if sum(dat2use)
    fname = datFiles(dat2use).name;
else
    error('no .dat file matching script search')
end

%let's duplicate the .dat file so you can keep the original
disp('Creating PC_removed.dat file')
copystring = ['! rsync -a -v -W -c -r ', ...
                        fname, ' ',...
                        [fname(1:end-4),'_PC_removed.dat']];
eval(copystring)
%update to do our changes on the new copied data
fname = [fname(1:end-4),'_PC_removed.dat'];


%starting script
disp('Removing 1st PC from data')
try
    %loading xml file
    d   = dir('*.xml');
    if ~isempty(d)
        par = LoadXml(fullfile(basepath,d(1).name));
    else
        error('the .xml file is missing')
    end
    %getting number of channels
    nbChan = par.nChannels;
    
    infoFile = dir(fname);
    
    chunk = 1e6;
    nbChunks = floor(infoFile.bytes/(nbChan*chunk*2));
    warning off
    if nbChunks==0
        chunk = infoFile.bytes/(nbChan*2);
    end
    
    for ix=0:nbChunks-1
        m = memmapfile(fname,'Format','int16','Offset',ix*chunk*nbChan*2,'Repeat',chunk*nbChan,'writable',true);
        d = m.Data;
        d = double(reshape(d,[nbChan chunk])');
        
        %doing for each shank separately with channels that were not skipped
        %on neuroscope
        new_d = zeros(size(d));
        for a = 1:length(par.AnatGrps)
            %channel number of this shank
            CHnum = par.AnatGrps(a).Channels+1; %neuroscope is 0-based
            %channels to use
            useCH=~logical(par.AnatGrps(a).Skip);
            
            if sum(useCH)>0 && sum(useCH)==length(useCH)
                
                [U,S,V] = svd(normalize(d(:,CHnum(useCH)),'zscore'),'econ');
                %reconstruct data
                new_d(:,CHnum(useCH)) = U(:,2:end)*S(2:end,2:end)*V(:,2:end)';
            end
        end
        new_d = new_d';
        m.Data = int16(new_d(:));
        clear d m
    end
    %close(h)
    
    
    newchunk = infoFile.bytes/(2*nbChan)-nbChunks*chunk;
    
    if newchunk
        m = memmapfile(fname,'Format','int16','Offset',nbChunks*chunk*nbChan*2,'Repeat',newchunk*nbChan,'writable',true);
        d = m.Data;
        d = double(reshape(d,[nbChan newchunk])');
        
        %doing for each shank separately with channels that were not skipped
        %on neuroscope
        new_d = zeros(size(d));
        for a = 1:length(par.AnatGrps)
            %channel number of this shank
            CHnum = par.AnatGrps(a).Channels+1; %neuroscope is 0-based
            %channels to use
            useCH=~logical(par.AnatGrps(a).Skip);
            
            if sum(useCH)>0 && sum(useCH)==length(useCH)
                
                [U,S,V] = svd(normalize(d(:,CHnum(useCH)),'zscore'),'econ');
                %reconstruct data
                new_d(:,CHnum(useCH)) = U(:,2:end)*S(2:end,2:end)*V(:,2:end)';
            end
        end
        new_d = new_d';
        m.Data = int16(new_d(:));
        clear d m
    end
    warning on
    returnVar = 1;
    msg = '';
    
catch
    fprintf(['Error occurred in processing ' fname '. File not processed.\n']);
    keyboard
    returnVar = 0;
    msg = lasterr;
end
clear m
end