function imgw = wedgify(img, r, theta, p, maskVal, addFix)

% imgw = wedgify(img, r, theta, p, maskVal, addFix)
%
% assumes only two phases of the cycle
% img in range [0 1]
% maskVal should be in img range

% maskVal = 0.5;

for iIm = 1:length(img)
    
    for iPhase = 1:2
        
        w = img{iIm};
        
        idx1 = r>p.wedgeRadial(1);
        idx2 = r<p.wedgeRadial(2);
        
        if rem(iPhase,2)==1 % cycle 1, 3, etc.
            
            switch p.mappingType
                
                case 'meridian'
                    % meridian mapping
                    idx3 = theta<p.wedgePolar/2;
                    idx4 = theta>2*pi-p.wedgePolar/2;
                    idx5 = theta<pi+p.wedgePolar/2;
                    idx6 = theta>pi-p.wedgePolar/2;
                    idx =  idx1.*idx2.*(idx3+idx4+(idx5.*idx6));
                    
                case'hemifieldUD'
                    % maps upper and lower vis field
                    idx5 = theta<3*pi/2+p.wedgePolar/2;
                    idx6 = theta>3*pi/2-p.wedgePolar/2;
                    idx =  idx1.*idx2*((idx5.*idx6));
                    
                case 'hemifieldLR'
                    % maps left and right hemifield
                    idx5 = theta<pi+p.wedgePolar/2;
                    idx6 = theta>pi-p.wedgePolar/2;
                    idx =  idx1.*idx2.*((idx5.*idx6));
                    
                otherwise
                    print 'Retinotopic mapping type not found.'
            end
            
            
        else % iPhase 0, 2, etc.
            
            switch p.mappingType
                
                case 'meridian'
                    % meridian mapping
                    idx3 = theta>pi/2-p.wedgePolar/2;
                    idx4 = theta<pi/2+p.wedgePolar/2;
                    idx5 = theta<3*pi/2+p.wedgePolar/2;
                    idx6 = theta>3*pi/2-p.wedgePolar/2;
                    idx =  idx1.*idx2.*((idx3.*idx4)+(idx5.*idx6));
                    
                case 'hemifieldUD'
                    % maps upper and lower vis field
                    idx3 = theta>pi/2-p.wedgePolar/2;
                    idx4 = theta<pi/2+p.wedgePolar/2;
                    idx =  idx1.*idx2.*((idx3.*idx4));
                    
                case 'hemifieldLR'
                    % maps left and right hemifield
                    idx3 = theta<p.wedgePolar/2;
                    idx4 = theta>2*pi-p.wedgePolar/2;
                    idx =  idx1.*idx2.*(idx3+idx4);
                    
                otherwise
                    print 'Retinotopic mapping type not found.'
            end
            
        end
        
        w(idx==0) = maskVal;
        
        if addFix
            w(r<p.fixSize) = 1; % add fixation
        end
        
        imgw{iIm, iPhase} = w;
        
    end
    
end