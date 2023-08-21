/////////////////////////////////////////////////////////////////// PROFILE SECTION => REFRESH BUTTON
/////////////////////////////////////////////////////////////////////////////////////////////////////
document.getElementById("profile-content__refresh-btn").onclick = function() {

    fetch("/refresh")
    .then(response => response.json())
    .then(data => {

        const profile_table = document.getElementById("profile-table");
        const profile_message = document.getElementById("profile-message");

        if (data.length == 0){
            profile_table.style.display = "none";
            profile_message.style.display = "block";
            return;
        }
        else {
            profile_table.style.display = "block";
            profile_message.style.display = "none";
        }

        let table = document.getElementById("profile-reports");
        table.innerHTML = "";

        data.forEach(problem => {

            let tr = document.createElement("tr");

            // NAME
            let td1 = document.createElement("td");
            td1.innerHTML = problem["name"];
            tr.append(td1);

            // METHOD
            let td2 = document.createElement("td");
            td2.innerHTML = problem["method"];
            tr.append(td2);

            // QUADRATURE
            let td3 = document.createElement("td");
            td3.innerHTML = problem["quadrature"];
            tr.append(td3);

            // MESH
            let td4 = document.createElement("td");
            td4.innerHTML = problem["mesh"];
            tr.append(td4);

            // TOLERANCE
            let td5 = document.createElement("td");
            td5.innerHTML = parseFloat(problem["tolerance"]).toExponential(2);
            tr.append(td5);

            // DATE
            let td6 = document.createElement("td");
            td6.innerHTML = problem["timestamp"];
            tr.append(td6);

            // LOAD BUTTON
            let td7 = document.createElement("td");
            let load = document.createElement("button");
            load.setAttribute("data-id", problem["id"]);
            load.className = "btn btn-success btn-sm";
            load.innerHTML = "<i class='bx bx-upload'></i>";
            load.onclick = function () {
                load_problem(this.dataset.id)
            }
            td7.append(load);
            tr.append(td7);

            // REMOVE BUTTON
            let td8 = document.createElement("td");
            let remove = document.createElement("button");
            remove.setAttribute("data-id", problem["id"]);
            remove.className = "btn btn-danger btn-sm";
            remove.innerHTML = "<i class='bx bx-trash' ></i>";
            remove.onclick = function () {
                remove_problem(this.dataset.id)
            }
            td8.append(remove);
            tr.append(td8);

            table.append(tr);

        })
    })
}


// LOAD PROBLEM
function load_problem(id){
    fetch(`load/${id}`)
    .then(response => response.json())
    .then(data => {
        if (data["returnCode"] == "ERROR"){
            console.log("Error loading problem.");
        }
        else {
  
            ZN = data["ZN"];
            ZON = data["ZON"];
            XREG = data["XREG"];
            XDOM = data["XDOM"];
            YREG = data["YREG"];
            YDOM = data["YDOM"];
            ZMAP = data["ZMAP"];
            QMAP = data["QMAP"];
            BC = data["BC"];
            QUAD = data["QUAD"]; 
            TOL = data["TOL"];
            METH = data["METH"];
            ITER = data["ITER"];
            CPU = data["CPU"];
            MFLUX = data["MFLUX"];
            MFLOW = data["MFLOW"];
            XFLOW = data["XFLOW"];
            YFLOW = data["YFLOW"];

            // REDIRECT TO RESULT SECTION
            document.getElementById("results-link").click();

            // PLOT RESULTS
            results_plot("plot_flux");

            // DISPLAY RESULTS IN CONSOLE
            console_process_cmd("clear");
            console_process_cmd("results");
            console_process_cmd("calculation_done");
            
        }
      
    });
}



// REMOVE PROBLEM
function remove_problem(id){

    // SEND DATA
    fetch(`remove/${id}`)
    .then(response => response.json())
    .then(data => {
  
      if (data["returnCode"] == "OK"){
  
        // REFRESH DATA
        document.getElementById("profile-content__refresh-btn").click();
  
      }
      
    });
  }

