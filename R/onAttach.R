
#[export special]
.onAttach <- function(lib, pkg) {
	version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),"Version")
	packageStartupMessage(paste("\n",pkg["name"],": ",version,sep = ""))
	msg <-  paste(
			r"( ___ __ __ __    __ __ __ __ _        _            __ __ __ __     __ __ __ __ __  )",
			r"(|  __ __ __  |  |  __ __ __ _/       / \          |  __ __ __ /   /__ __ __ __ __\ )",
			r"(| |        | |  | |                 / _ \         | |                   / /        )",
			r"(| |        | |  | |                / / \ \        | |                  / /         )",
			r"(| |__ __ __| |  | |__ __ __       / /   \ \       | |__ __ __ _       / /_/\       )",
			r"(|    __ __ __|  |  __ __ __|     / /__ __\ \      |_ __ __ _   |     / __  /       )",
			r"(|   \           | |             / _ _ _ _ _ \                | |     \/ / /        )",
			r"(| |\ \          | |            / /         \ \               | |       / /         )",
			r"(| | \ \__ _ _   | |           / /           \ \     _ __ __ _| |      / /          )",
			r"(|_|  \__ __ _\  |_|          /_/             \_\   /_ __ __ ___|      \/         team)"
	,sep="\n")
	packageStartupMessage(msg)
}