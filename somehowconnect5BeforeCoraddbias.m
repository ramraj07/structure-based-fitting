if 1==1
    iso = readtable('final_original.csv');
    isoo = iso;
    %%
     iso = isoo( isoo.Var4<1000 & isoo.Var4>-1000 ...
         & isoo.Var3<10000 & isoo.Var3>-10000,:);
  
    scatter3(iso.Var2,iso.Var3,iso.Var4,1)
    associations = iso.Var1*0;
    startpoint = [8484 14395 12389 10202];
    %  startpoint = [4619];[ 5323 1272 ];
    thispoints = startpoint;
    xvals = iso.Var2;
    yvals = iso.Var3;
    zvals = iso.Var4;
    
    %%
    hold off;
    clear ncomps
    nc = 300;
    rng('default');
    [idx,c] = kmeans([xvals yvals zvals],nc);
    cmap = hsv(nc);
    rr = randperm( max(idx(:)));
    for ii=1:numel(rr)
        i = rr(ii);
        scatter3(xvals(idx==i),...
            yvals(idx==i),...
            zvals(idx==i),1,cmap(ii,:));
        hold on;
        ncomps(i) = sum(idx==i);
    end
    corig = c;
    view(-10,90);
    %%
%      cla
%     hold off
%     c = corig(ncomps>40,:);
% %     scatter3(c(:,1),c(:,2),c(:,3),ncomps,associations)
%     scatter3(c(:,1),c(:,2),c(:,3),ncomps(ncomps>40))
%     for i=1:size(c,1)
%         text(c(i,1)+10,c(i,2)+10,c(i,3),num2str(i));
%     
%     end
    %%
     tubules = {...
           [201 223 61 132 289 88 88 69 47 246 261 103 53 208 115 293 117 155 290 82 25 183 113 279 12 30 ...
           236 30 242 281 284 59 174 217 142 210 277 87 44 131 123 153 51 33 271 189 105 17 64 197 179 39 ...
           296 120 188 90 163 255 298 291 199 4 157 127 229 57 278 211 164 18 234 219 299 151 38 73 167 104 226 15 254 250 243 161 185 130 176 285 94 31 263 258 21],...
        [32 85 215 251 77 231 7 116 84 160 128 96 170 205 273 158 206 37 48 178 260 114 180 3 214 137 177 80 ...
        9 268 280 264 232 22 227 20 207 78 262 288 1 203 235 52 213 272 209 106 19 238 186 101 184 26 240 83  ...
        72 252 224 97 46 150 187 74 152 168 126 247 259 35 191 67 111 45 166 300 108 24 286 245 228 10 100 266 ...
        148 40 283 239 60 68 70 175 190 13 182 86 241 76 93 28 216 71 282 122 267 36 193 109 256 135 134 125 ...
        241 237],...      
        [221 63 248 8 133 43 233 143 56 102 29 89 212 162 23 62 49 138 192 149 204 222 292 156 75 269 2 119 ...
        99 172 129 65 91 34 145 112 154 171 287 165 173 202 275 140 198 141 41 244 136 118 121 181 5 139 270 ...
        107 58 146 79 50 295 249 124 220 14 274 225 257 95 110 194 66 230 276 27 218 16 54 196 200 297 147 294 ...
        98 144 55 81 159 92 11 195 42 265 253 169 6]};
    tubules{1} = [tubules{1}   ];
    tubules{2} = [tubules{2} ];
    tubules{3} = [tubules{3} ];
%     tubules = {[],[],[]};
    colors = zeros(size(corig,1),1);
    for j=1:numel(tubules)
        colors(tubules{j}) = j;
    end
    cla
    hold off
    scatter3(corig(:,1),corig(:,2),corig(:,3),(ncomps),colors)
    for j=1:size(corig,1)
        if colors(j)==0
             text(corig(j,1),corig(j,2),corig(j,3),num2str(j));
        end
    end
    
    % scatter3(c(:,1),c(:,2),c(:,3),ncomps,associations)
    
    %%
    belongstoatubule = ismember(idx,cell2mat(tubules));
    clustertubuleids = zeros(size(c,1),1);
    for j=1:numel(tubules)
        tubuleidentitymat(:,j) = ismember(idx,tubules{j});
        clustertubuleids(ismember(1:size(c,1),tubules{j})) = j;
    end
    for j=1:size(tubuleidentitymat,1)
        this = find(tubuleidentitymat(j,:),1,'first');
        if isempty(this)
            tubuleidentity(j) = 0;
        else
            tubuleidentity(j) = this;
        end
    end
    tubuleid = xvals*0;
    for i=1:numel(xvals)
        if belongstoatubule(i)
            tubuleid(i) = tubuleidentity(i);
        else
            
            tubuleid(i) = 0;
        end
        
    end
    scatter3(xvals,yvals,zvals,3,tubuleid)
end
%%
%%
hold off
zvalscorrected = zvals;
xvalscorrected = xvals;
yvalscorrected = yvals;
zvalsfitted = zvals;
xvalsfitted = xvals;
yvalsfitted = yvals;
for j=1:numel(tubules)
    %%
    xb = xvals(tubuleid==j);
    yb = yvals(tubuleid==j);
    zb = zvals(tubuleid==j)-12.76;
    
    nc = 200;
    [idx2{j},clustercenters{j}] = kmeans([xb yb ],nc);
    cmap = hsv(nc);
    hold off
    rr = randperm( max(idx2{j}(:)));
    zbcorrected =zb;
    xbcorrected = xb;
    ybcorrected = yb;
    for ii=1:numel(rr)
        i = rr(ii);
        zbthis = zb(idx2{j}==i);
        zbthis = (zbthis*0)+median(zbthis);
        zbcorrected(idx2{j}==i) = zbthis;
        xbcorrected(idx2{j}==i) = 0*zbthis+median(xb(idx2{j}==i));
        ybcorrected(idx2{j}==i) = 0*zbthis+median(yb(idx2{j}==i));
  
        ncomps(i) = sum(idx2{j}==i);
    end
    % hold off;
    % scatter3(xb,yb,zb,1)
    % hold on;
    %     scatter3(xb,yb,zbcorrected,1)
    %
    %     view(0,0);
    zvalscorrected(tubuleid==j) = zbcorrected;
    yvalscorrected(tubuleid==j) = ybcorrected;
    xvalscorrected(tubuleid==j) = xbcorrected;
    %     distmat = squareform(pdist([xb yb zb]));
    %     [~,farthestpt] = max(max(distmat));
    %     orderofappearance = xb*0;
    %     currentpoint = farthestpt;
    %     for i=1:numel(xb)
    %         orderofappearance(currentpoint) = max(orderofappearance)+1;
    %         thesedistances = distmat(:,currentpoint);
    %         thesedistances(orderofappearance~=0) = 1e10;
    %         [jumpsize,currentpoint] = min(thesedistances);
    %         if jumpsize>1000
    %             break;
    %         end
    %
    %     end
    %     zf = zvals(orderofappearance~=0);
    %     plot(zf(orderofappearance(orderofappearance~=0)))
    %     [~,sortorder] = sort(orderofappearance);
    %      scatter3(xvals(tubuleid==j),yvals(tubuleid==j),zvals(tubuleid==j),1,orderofappearance);
    fitobject = fit([xb,yb],zbcorrected,'lowess','span',0.30);
    zbfitted =fitobject(xb,yb);
    zvalsfitted(tubuleid==j) = zbfitted;
    
%     fitobjectx = fit([yb,zb],xbcorrected,'lowess','span',0.30);
%     xbfitted =fitobjectx(yb,zb);
%     xvalsfitted(tubuleid==j) = xbfitted;
%     fitobjecty = fit([xb,zb],ybcorrected,'lowess','span',0.30);
%     ybfitted =fitobjecty(xb,zb);
%     yvalsfitted(tubuleid==j) = ybfitted;
end
zvalsfittedorig = zvalsfitted;
xvalscorrectedorig = xvalscorrected;
yvalscorrectedorig = yvalscorrected;
%%
zvalsfitted = zvalsfittedorig;
xvalscorrected = xvalscorrectedorig;
yvalscorrected = yvalscorrectedorig;
 switchback = abs(xvalscorrected-xvals)+abs(yvalscorrected-yvals)<40;
 xvalscorrected(switchback)=xvals(switchback);
 yvalscorrected(switchback)=yvals(switchback);
if 1==0
    %%
%     close all
    
%     scatter3(xb,yb,zb);

end
% for j=1:numel(zvalsfitted)
%     [~,minpos] = min(abs(lookupvals-zvalsfitted(j)));
% %    zvalsfitted(j) = zvalsfitted(j)-zcorrection(minpos);
% end
% %%
% overlapx = ove.Var2;
% overlapy = ove.Var3;
% overlapz = ove.Var4;
% correctedoverlapz = overlapz;
% for i=1:numel(overlapz)
%     distances = abs(c(:,1)-overlapx(i))+...
%         abs(c(:,2)-overlapy(i));%+        abs(c(:,3)-overlapz(i));
%     [~,closestOriginalCluster] = min(distances);
%     correctedoverlapz(i) = median(zvalscorrected(idx==closestOriginalCluster));
% %     tubuleidentity = clustertubuleids(closestOriginalCluster);
% %     if tubuleidentity==0
% %         continue;
% %     end
% %     distances2 =  abs(clustercenters{tubuleidentity}(:,1)-overlapx(i))+...
% %         abs(clustercenters{tubuleidentity}(:,2)-overlapy(i));
% %     [~,closestTubuleCluster] = min(distances2);
% %     correctedoverlapz(i) = zcorrectedperclusterpertubule{tubuleidentity}(closestTubuleCluster);
% end
%
 close all

% % scatter3(xvals,yvals,zvals,3);
 scatter3(gt.xnano,gt.ynano,gt.znano,1);
 hold on;
%   scatter3(xvals,yvals,zvals,1);
 ax = axis;
  hold on;
 scatter3(xvalscorrected,yvalscorrected,zvalsfitted,1);
 hold on;
%   scatter3(xvals,yvals,zvals,1);
% scatter3(ove.Var2,ove.Var3,correctedoverlapz,1);
% view(-90,90)
axis(ax);
%
isoc = iso;
drawnow;
%%
% 
% ovec = ove;
% 
% ovec2 = isoc(1:size(ovec,1),:);
% ovec2 = ovec;
% combined = [isoc; ovec2];
% ovec2.Var4 = correctedoverlapz;
isoc.Var4 = zvalscorrected;
% combinedc = [isoc; ovec2];
writetable(iso,'combined-uncorrected-z.csv');
writetable(isoc,'combined-corrected-z.csv');
isoc.Var4 = zvalsfitted;
isoc.Var2 = xvalscorrected;
isoc.Var3 = yvalscorrected;
writetable(isoc,'combined-corrected-fitted-z.csv');
disp('written')
%%
save(['runat',datestr(now,30)]);
