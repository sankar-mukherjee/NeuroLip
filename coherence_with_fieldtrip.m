
clear;clc;

data_path = '\\pcespvicon\projects\current\listen_italian_motor_entrainment\analysis\python/data/coherence_analysis_matlab/';
freq_band = [0.5:0.5:15];
subject_name = {'Alice','Andrea','Daniel','Elena','Elenora','Elisa','Federica','Francesca','Gianluca1','Giada','Giorgia',...
    'Jonluca','Laura','Leonardo','Linda','Lucrezia','Manu','Marco','Martina','Pagani','Pasquale','Sara',...
    'Silvia','Silvia2','Tommaso'};

feature = {'envelop';'lipaparature'};
delay = [-0.5:0.025:0.5];
frsmooth=2;
triallen=2;
method = 'coh';
pad = 6; %'nextpow2'
XXX = [];
XXX_S = [];
surrogate=0;
fsample = 100;

for d=1:length(delay)
    fd=[];
    
    for sub=1:length(subject_name)
        if(d==21)
            load([data_path subject_name{sub} '-trialLength' num2str(triallen) '-delay0.0.mat']);
        else
            load([data_path subject_name{sub} '-trialLength' num2str(triallen) '-delay' num2str(delay(d)) '.mat']);
        end
        
        % change label
        A =[];
        for i=1:59
            A{i} = strtrim(label(i,:));
        end
        label=A';
        label{60} = 'envelop';
        label{61} = 'lipaparature';
        
        %% data
        % fieldtrip data
        X = [];
        X.fsample = fsample;
        X.label= label;
        for t=1:size(data,1)
            X.trial{t} = squeeze(data(t,:,:));
            X.time{t} = [0:1/X.fsample:triallen];            
        end
        
        
        AA = [];
        for ff=1:2
            cfg            = [];
            cfg.output     = 'powandcsd';
            cfg.method     = 'mtmfft';
            cfg.taper      =  'dpss' ;
            cfg.foi    = freq_band;
            cfg.channelcmb = {'EEG' feature{ff}};
            cfg.pad        = pad;
            cfg.tapsmofrq  = frsmooth;
            freqA           = ft_freqanalysis(cfg, X);

            cfgg            = [];
            cfgg.method     = method;
            % cfgg.complex     = metric_value;
            a = ft_connectivityanalysis(cfgg, freqA);
            if(method=='coh')
                AA{ff,1} = a.cohspctrm;AA{ff,2} = round(a.freq,2);   
            else
                AA{ff,1} = a.plvspctrm;AA{ff,2} = round(a.freq,2);   
            end
            
            
            freq = round(a.freq,1);            
        end
        
        %% surrogate
        if(surrogate)
            XX = [];
            for chh=1:59
                a = randperm(size(data,1));
                data_S = cat(2,data(:,1:59,:),data(a,60:end,:));
                % fieldtrip data

                for t=1:size(data,1)
                    X.trial{t} = squeeze(data_S(t,:,:));                
                end


                for ff=1:2
                    cfg            = [];
                    cfg.output     = 'powandcsd';
                    cfg.method     = 'mtmfft';
                    cfg.taper      =  'dpss' ;
                    cfg.foi    = freq_band;
                    cfg.channelcmb = {'EEG' feature{ff}};
                    cfg.pad        = pad;
                    cfg.tapsmofrq  = frsmooth;
                    freqA           = ft_freqanalysis(cfg, X);

                    cfgg            = [];
                    cfgg.method     = 'coh';
                    % cfgg.complex     = metric_value;
                    a = ft_connectivityanalysis(cfgg, freqA);
                    XX{chh,ff} = a.cohspctrm(chh,:);
                end
            end
            AA{1,2} = cell2mat(XX(:,1));        AA{2,2} = cell2mat(XX(:,2));
            fd{sub,3}   = AA{1,2};
            fd{sub,4}   = AA{2,2};
        end
       
        fd{sub,1}   = AA{1,1};
        fd{sub,2}   = AA{2,1};
        
        fd{sub,3}   = AA{1,2};
		fd{sub,4}   = AA{2,2};
    end
    
    XXX{d,1} = cat(3,fd{:,1});
	XXX{d,2} = cat(3,fd{:,2});
	
	XXX{d,3} = cat(1,fd{:,3});
	XXX{d,4} = cat(1,fd{:,4});
    
    if(surrogate)
        XXX_S{d,1} = cat(3,fd{:,3});XXX_S{d,2} = cat(3,fd{:,4});
    end
    clc;
    disp(d);
end

A = cat(4,XXX{:,1});
B = cat(4,XXX{:,2});

envelop = permute(A,[3 1 2 4]);
lip = permute(B,[3 1 2 4]);


AA = cat(4,XXX{:,3});
BB = cat(4,XXX{:,4});

envelop_freq = squeeze(permute(AA,[3 1 2 4]));
lip_freq = squeeze(permute(BB,[3 1 2 4]));

if(surrogate)
    A = cat(4,XXX_S{:,1});
    B = cat(4,XXX_S{:,2});

    envelop_S = permute(A,[3 1 2 4]);
    lip_S = permute(B,[3 1 2 4]);
else
    lip_S=0;
    envelop_S=0;
end

time = [0:1/X.fsample:triallen];
% load([data_path subject_name{sub} '-delay' num2str(delay(d)) '.mat']);

save(['\\PCESPVICON\projects\current\listen_italian_motor_entrainment\analysis\python\data\' method '_fieldtrip_pad_' ...
    num2str(pad) '_' num2str(frsmooth) 'hzsmoothDPSS_' num2str(freq_band(1))...
    '-' num2str(freq_band(end)) 'Hz-' num2str(triallen) '.mat'],'envelop','lip','envelop_S',...
	'envelop_freq','lip_freq','lip_S','freq','time','label','delay');



