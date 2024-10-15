%[dont read]
% System Rd macros

% These macros are automatically loaded whenever R processes an Rd file.

% Packages may define their own macros, which are stored in man/macros/*.Rd in
% the source, help/macros/*.Rd after installation.  Those will be processed after
% this file but before every Rd file in a package.

% Packages may request inclusion of macros from other packages using the 
% LoadRdMacros line in the DESCRIPTION file, e.g.
%  LoadRdMacros:  pkgA 
% These are loaded after the system macros and before the current package macros.

% Individual Rd files may define their own macros.



% To get the package author at build time from the DESCRIPTION file
\newcommand{\packagePackage}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1")[["Package"]]}}
\newcommand{\packageType}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1")[["Type"]]}}
\newcommand{\packageVersion}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1")[["Version"]]}}
\newcommand{\packageDate}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1")[["Date"]]}}
\newcommand{\packageLicense}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1")[["License"]]}}
\newcommand{\packageName}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1")[["Package"]]}}
\newcommand{\packageAuthors}{\Sexpr[results=rd,stage=build]{utils::packageDescription("#1")[["Author"]]}}