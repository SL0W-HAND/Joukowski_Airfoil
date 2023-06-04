%% =============================================================================
%% ppgcaposter.cls
%%
%% Copyright 2017 Adriano J. Holanda, DCM/FFCLRP/USP
%%
%% This software is distributed under MIT License
%% https://opensource.org/licenses/MIT
%% =============================================================================

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{Hylangtechposter}[2017/11/27 class to DCM/FFCLRP/USP graduate manuscript template]

% --- Class structure: declaration of options part
% ---
% This class extends the a0poster class
% Read all the documentclass options; pass them to a0poster,
% unless the file "<currentoption>.ppgcaposter" exists, then it is loaded
\DeclareOption*{\InputIfFileExists{\CurrentOption.Hylangtechposter}{}{
    % 
    \PassOptionsToClass{\CurrentOption}{a0poster}}}

% --- Class structure: execution of options part
% ---
\ProcessOptions \relax
% --- Class structure: declaration of options part
% ---
\LoadClass{a0poster}

\ExecuteOptions{final}
\ProcessOptions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PACKAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{fontspec}
\RequirePackage{multicol} % This is so we can have multiple columns of text side-by-side
\columnsep=50pt % This is the amount of white space between the columns in the poster
\columnseprule=3pt % This is the thickness of the black line between the columns in the poster
\usepackage[spanish]{babel}

\RequirePackage[svgnames]{xcolor}
\RequirePackage{amsmath}
\RequirePackage{graphicx} % Required for including images
\RequirePackage{booktabs} % Top and bottom rules for table
\RequirePackage[font=small,labelfont=bf]{caption} % Required for specifying captions to tables and figures

\RequirePackage{amsfonts, amsmath, amsthm, amssymb} % For math fonts, symbols and environments
\RequirePackage{wrapfig} % Allows wrapping text around tables and figures
\RequirePackage{tikz}
\usetikzlibrary{calc}
\usetikzlibrary{fadings}
\RequirePackage{caption}
\RequirePackage{etoolbox}
\DeclareGraphicsRule{.ai}{pdf}{.ai}{}
\RequirePackage{lipsum}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FONTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% size 24.88pt
\renewcommand{\tiny}{\fontsize{12}{14}\selectfont}
\renewcommand{\scriptsize}{\fontsize{14.4}{18}\selectfont}
\renewcommand{\footnotesize}{\fontsize{17.28}{22}\selectfont}
\renewcommand{\small}{\fontsize{20.74}{25}\selectfont}
\renewcommand{\normalsize}{\fontsize{24.88}{30}\selectfont}
\renewcommand{\large}{\fontsize{29.86}{37}\selectfont}
\renewcommand{\Large}{\fontsize{35.83}{45}\selectfont}
\renewcommand{\LARGE}{\fontsize{43}{54}\selectfont}
\renewcommand{\huge}{\fontsize{51.6}{64}\selectfont}
\renewcommand{\Huge}{\fontsize{61.92}{77}\selectfont}
\renewcommand{\veryHuge}{\fontsize{74.3}{93}\selectfont}
\renewcommand{\VeryHuge}{\fontsize{89.16}{112}\selectfont}
\renewcommand{\VERYHuge}{\fontsize{107}{134}\selectfont}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A0 dimensions 841 x 1189 mm  33.1 x 46.8 in
%\special{papersize=46.8in,33.1in}
\setlength{\paperwidth}{92cm}
\setlength{\paperheight}{169cm}
\setlength{\textwidth}{85cm}
\setlength{\textheight}{163cm}

\topskip0pt
\setlength{\headheight}{0pt}
\setlength{\headsep}{0pt}
\setlength{\topmargin}{0pt}
\setlength{\oddsidemargin}{0pt}

\setlength{\oddsidemargin}{1cm}

%\setlength{\columnsep}{2cm}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMANDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand\printheader{
  \begin{minipage}[ht]{0.15\linewidth}
    \includegraphics[scale=1]{uach logo.png}
  \end{minipage}\vspace{2cm}
  \begin{minipage}[ht]{0.45\linewidth}\bf
    {\huge Universidad Autónoma de Chihuahua}\\[0.4cm]
    {\huge Facultad de Ingenería}\\[0.3cm]
    \Large Ingenería Física
  \end{minipage} \hfill
  \begin{minipage}[ht]{0.21\linewidth}\bf
    \centering    \includegraphics[scale=2.1]{logo ingenieria.png} \\[-1cm]

  \end{minipage}
  %\begin{minipage}[ht]{0.2\linewidth}\bf
  %    \includegraphics[trim={0mm 5mm 5mm 2mm}, clip, scale=.35]{eu.pdf} \quad 
  %    \includegraphics[trim={32mm 55mm 30mm 5mm}, clip, scale=.63]{erc.ai} 
  %\end{minipage}
}