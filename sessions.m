Sessions.session = [];
Sessions.subject = {};

% summer sessions
% Sessions.session = [repmat(datetime(2022,7,15),1,2) datetime(2022,8,20)];
% Sessions.subject = {'237','238','269'};

% 2022 Dec cohort
% Sessions.session = [Sessions.session repmat(datetime(2022,12,4),1,3)];
% Sessions.subject = [Sessions.subject repmat({'001','002','277'},1,1)];


% 2023 Feb cohort
% Sessions.session = [Sessions.session repmat(datetime(2023,2,14),1,4)];
% Sessions.subject = [Sessions.subject repmat({'273','274','279','280'},1,1)];
% Sessions.session = [Sessions.session repmat(datetime(2023,2,14),1,2)];
% Sessions.subject = [Sessions.subject repmat({'279','280'},1,1)];


% 2023 Apr cohort
% Sessions.session = [Sessions.session repmat(datetime(2023,4,19),1,4)];
% Sessions.subject = [Sessions.subject repmat({'302','304','305','306'},1,1)];

% 2024 ablation mice
Sessions.session = [Sessions.session repmat(datetime(2024,4,18),1,3) datetime(2024,7,5)];
Sessions.subject = [Sessions.subject repmat({'461','462','463'},1,1) '913'];

% 2024 control mice
Sessions.session = [Sessions.session repmat(datetime(2024,7,5),1,4)];
Sessions.subject = [Sessions.subject repmat({'903','905','906','907'},1,1)];

if ~exist('sessions_oi','var')
    sessions_oi = 1:numel(Sessions.session);
end