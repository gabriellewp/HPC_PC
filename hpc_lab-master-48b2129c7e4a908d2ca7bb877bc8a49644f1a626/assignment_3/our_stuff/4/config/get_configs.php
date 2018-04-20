<?php
	for($d1 = 8; $d1<=8192; $d1++){
		$mesh = 1/($d1-1);
		$d1_used = $d1-2;
		echo "<table><tr>";
		echo "<th style='width:250px'>".$mesh."</th>";
		echo "<th style='width:80px'>".$d1."</th>";
		for($colPerP = 1; $colPerP<=$d1_used; $colPerP++){
			$conf = ($d1_used/$colPerP)+1;
			if($conf == (int)$conf){
				// valid conf
				echo "<td style='width:40px;'>".$conf."</td>";
			}
		}
		echo "</tr></table>";
	}
?>