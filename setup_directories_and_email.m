function [] = setup_directories_and_email(platform)
% Setup parameters to email results.
%
% Args:
%   platform: String, either 'mac' or 'linux'.
%
% Returns:
%   none

% Import GPstuff and set paths.
if strcmp(platform, 'mac');  % Mac version.
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/GPstuff-4.6/
    matlab_install
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/Programs/
    run_mex_commands();
    cd ~/Google' Drive'/0-LIZHEN' RESEARCH'/gp/
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Programs')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Functions')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/Smoothing')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/convex-function')
    addpath('~/Google Drive/0-LIZHEN RESEARCH/gp/mbcr')
elseif strcmp(platform, 'linux');  % Linux version.
    cd ~/Documents/gp/GPstuff-4.6/
    matlab_install
    cd ~/Documents/gp/Programs/
    run_mex_commands();
    cd ~/Documents/gp/
    addpath('~/Documents/gp/')
    addpath('~/Documents/gp/Programs')
    addpath('~/Documents/gp/Functions')
    addpath('~/Documents/gp/Smoothing')
    addpath('~/Documents/gp/convex-function')
    addpath('~/Documents/gp/mbcr')
end

% Setup email parameters.
myaddress = 'eltegog@gmail.com';
myp = 'T0g.eltegog';
setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',myp);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

end

