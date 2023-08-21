////////////////////////////////////////////////////////////////////////// SETUP SECTION => VARIABLES
/////////////////////////////////////////////////////////////////////////////////////////////////////
let ZN = 0;
let ZON = [];
let XREG = 0;
let XDOM = [];
let YREG = 0;
let YDOM = [];
let ZMAP;
let QMAP;
let BC = {
    "left": 0,
    "right": 0,
    "top": 0,
    "bottom": 0
}



///////////////////////////////////////////////////////////////////// SETUP SECTION => INITIALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////
function setup_section_init() {

    // CLEAR UPLOAD GEOMETRY
    setup_upload_file_reset();

}



//////////////////////////////////////////////////////////////////////// SETUP SECTION => UPLOAD FILE
/////////////////////////////////////////////////////////////////////////////////////////////////////

// FUNCTION FOR LOAD EXTERNAL FILE
function setup_upload_file() {
    let input = document.getElementById("setup-section__input--upload-file");
    let files = input.files;
    if(files.length == 0) return;
    const file = files[0];
    let reader = new FileReader();
    reader.onload = function (e){
        let input = document.getElementById("setup-section__input--upload-file");
        let err = document.getElementById("setup-section__input--upload-file__msg");
        let data;
        try {
            data = JSON.parse(e.target.result);
        }
        catch(e) {
            input.classList.remove("is-valid");
            input.classList.remove("is-invalid");
            input.classList.add("is-invalid");
            err.className = "invalid-feedback";
            err.innerHTML = "Invalid file";
            return;
        }
        if (setup_upload_file_validator(data)){

            // LOAD GEOMETRY
            ZN = data["ZN"];
            ZON = data["ZON"];
            XREG = data["XREG"];
            XDOM = data["XDOM"];
            YREG = data["YREG"];
            YDOM = data["YDOM"];
            ZMAP = data["ZMAP"];
            QMAP = data["QMAP"];
            BC = data["BC"];

            // RESET SOLUTION AND RESULTS
            solution_quadrature_reset();
            results_reset();
            if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
            if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
            if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
            document.getElementById("initial_img").style.display = "block";

            input.classList.remove("is-valid");
            input.classList.remove("is-invalid");
            input.classList.add("is-valid");
            err.className = "valid-feedback";
            err.innerHTML = "Upload Completed";
        }
        else {
            input.classList.remove("is-valid");
            input.classList.remove("is-invalid");
            input.classList.add("is-invalid");
            err.className = "invalid-feedback";
            err.innerHTML = "Invalid data";
        }
    }
    reader.readAsText(file);
}

// FUNCTION FOR RESET FILE INPUT
function setup_upload_file_reset() {
    let input = document.getElementById("setup-section__input--upload-file");
    let err = document.getElementById("setup-section__input--upload-file__msg");
    input.value = "";
    input.classList.remove("is-valid");
    input.classList.remove("is-invalid");
    err.className = "";
    err.innerHTML = "";
}

// FUNCTION TO CHECK FOR INVALID DATAFILE
function setup_upload_file_validator(data){
    
    // ZONES
    if (data["ZN"] <= 0 || data["ZN"] == undefined) return false;
    if (data["ZON"] == undefined || data["ZON"].length != data["ZN"]) return false;
    for (let z = 0; z < data["ZON"].length; z++){
        if (data["ZON"][z].st <= 0 || data["ZON"][z].st == undefined) return false;
        if (data["ZON"][z].ss < 0 || data["ZON"][z].ss == undefined || data["ZON"][z].ss >= data["ZON"][z].st) return false;
    }

    // X REGIONS
    if (data["XREG"] <= 0 || data["XREG"] == undefined)  return false;
    if (data["XDOM"] <= undefined || data["XDOM"].length != data["XREG"]) return false;
    for (let i = 0; i < data["XREG"]; i++){
        if (data["XDOM"][i].len <= 0 || data["XDOM"][i].len == undefined)  return false;
        if (data["XDOM"][i].nc < 0 || data["XDOM"][i].nc == undefined || !Number.isInteger(data["XDOM"][i].nc))  return false;
    }

    // Y REGIONS
    if (data["YREG"] <= 0 || data["YREG"] == undefined)  return false;
    if (data["YDOM"] <= undefined || data["YDOM"].length != data["YREG"]) return false;
    for (let i = 0; i < data["YREG"]; i++){
        if (data["YDOM"][i].len <= 0 || data["YDOM"][i].len == undefined)  return false;
        if (data["YDOM"][i].nc < 0 || data["YDOM"][i].nc == undefined || !Number.isInteger(data["YDOM"][i].nc))  return false;
    }

    // ZMAP
    if (data["ZMAP"] == undefined || data["ZMAP"].length != data["YREG"])  return false;
    if (data["ZMAP"][0] == undefined || data["ZMAP"][0].length != data["XREG"])  return false;
    for (let i = 0; i < data["XREG"]; i++){
        for (let j = 0; j < data["YREG"]; j++){
            if (data["ZMAP"][j][i] < 0 || data["ZMAP"][j][i] >= data["ZN"] || !Number.isInteger(data["ZMAP"][j][i])) return false;
        }
    }
    
    // QMAP
    if (data["QMAP"] == undefined || data["QMAP"].length != data["YREG"])  return false;
    if (data["QMAP"][0] == undefined || data["QMAP"][0].length != data["XREG"])  return false;
    for (let i = 0; i < data["XREG"]; i++){
        for (let j = 0; j < data["YREG"]; j++){
        if (data["QMAP"][j][i] < 0) return false;
        }
    }

    // BC
    if (data["BC"] == undefined) return false;
    if (data["BC"].left == undefined || (data["BC"].left < 0 && data["BC"].left != -1)) return false;
    if (data["BC"].right == undefined || (data["BC"].right < 0 && data["BC"].right != -1)) return false;
    if (data["BC"].top == undefined || (data["BC"].top < 0 && data["BC"].top != -1)) return false;
    if (data["BC"].bottom == undefined || (data["BC"].bottom < 0 && data["BC"].bottom != -1)) return false; 

    return true;
}



///////////////////////////////////////////////////////////////////// SETUP SECTION => ZONE SETTINGS
////////////////////////////////////////////////////////////////////////////////////////////////////

// FUNCTION FOR INITIALIZE ANY AVAILABLE ZONE
function setup_zone_init(){
    document.querySelector("#modal-zone__table").innerHTML = "";
    for (let z = 0; z < ZN; z++){
       setup_zone_load(z, ZON[z].st, ZON[z].ss);
    }
}

// FUNCTION FOR ADDING A NEW ZONE
function setup_zone_add(){
    setup_zone_load(ZN,"","");
    document.querySelector(`#st${ZN}`).disabled = false;
    document.querySelector(`#ss${ZN}`).disabled = false;
    let btn = document.querySelector(`#save${ZN}`);
    btn.dataset.action = "save";
    btn.className = "btn btn-success btn-sm";
    btn.innerHTML = "<i class='bx bx-save' ></i>";
    ZN = ZN + 1;
    let z = {
      "st": "",
      "ss": ""
    }
    ZON.push(z);
}

// FUNCTION FOR LOAD ZONES MATERIAL PARAMETERS
function setup_zone_load(z_id, z_st, z_ss){
    let table = document.getElementById("modal-zone__table");
    let tr = document.createElement("tr");
    tr.id = `row${z_id}`;
  
    // ID
    let id = document.createElement("td");
    id.className = "input-group-sm";
  
    let id_input = document.createElement("input");
    id_input.id = `zone_id${z_id}`;
    id_input.disabled = true;
    id_input.className = "form-control";
    id_input.type = "text";
    id_input.value = `${z_id + 1}`;
    id.append(id_input);
    tr.append(id);
  
    // TOTAL CROSS SECTION
    let st = document.createElement("td");
    st.className = "input-group-sm";
    let st_input = document.createElement("input");
    let st_err = document.createElement("div");
    st_input.className = "form-control";
    st_input.id = `st${z_id}`;
    st_err.id = `sterr${z_id}`;
    st_input.disabled = true;
    st_input.type = "text";
    st_input.setAttribute("autocomplete", "off");
    st_input.value = z_st;
    st.append(st_input);
    st.append(st_err);
    tr.append(st);
  
    // SCATTERING CROSS SECTION
    let ss = document.createElement("td");
    ss.className = "input-group-sm";
    let ss_input = document.createElement("input");
    let ss_err = document.createElement("div");
    ss_input.className = "form-control";
    ss_input.id = `ss${z_id}`;
    ss_err.id = `sserr${z_id}`;
    ss_input.disabled = true;
    ss_input.type = "text";
    ss_input.setAttribute("autocomplete", "off");
    ss_input.value = z_ss;
    ss.append(ss_input);
    ss.append(ss_err);
    tr.append(ss);
  
    // EDIT AND SAVE BUTTONS
    let td = document.createElement("td");
    let btn = document.createElement("button");
    btn.id = `save${z_id}`;
    btn.className = "btn btn-secondary btn-sm";
    btn.setAttribute("data-id", `${z_id}`);
    btn.setAttribute("data-action", "edit");
    btn.innerHTML = "<i class='bx bx-edit'></i>";
    btn.onclick = function(){
        if (setup_zone_validator(this.dataset.id) === true){
            let id = parseInt(this.dataset.id);
            let st = document.getElementById(`st${id}`);
            let ss = document.getElementById(`ss${id}`);
            let action = this.dataset.action;
            if (action == "edit"){
                st.disabled = false;
                ss.disabled = false;
                this.dataset.action = "save";
                this.className = "btn btn-success btn-sm";
                this.innerHTML = "<i class='bx bx-save' ></i>";
            }
            else {
                ZON[id].st = parseFloat(st.value);
                ZON[id].ss = parseFloat(ss.value);
                st.disabled = true;
                ss.disabled = true;
                this.dataset.action = "edit";
                this.className = "btn btn-secondary btn-sm";
                this.innerHTML = "<i class='bx bx-edit'></i>";
            }

            // RESET ZMAP
            ZMAP = undefined;
            
            // RESET SOLUTION AND RESULTS
            solution_quadrature_reset();
            results_reset();
            if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
            if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
            if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
            document.getElementById("initial_img").style.display = "block";
        }
    };
    td.append(btn);
    tr.append(td);
  
    // REMOVE BUTTON
    let td2 = document.createElement("td");
    let btn2 = document.createElement("button");
    btn2.id = `remove${z_id}`;
    btn2.className = "btn btn-danger btn-sm";
    btn2.setAttribute("data-id", `${z_id}`);
    btn2.innerHTML = "<i class='bx bxs-trash'></i>";
    btn2.addEventListener("click", function(){
  
        let id = parseInt(this.dataset.id);
        document.getElementById(`row${id}`).remove();
        ZON.splice(id, 1);
        for(let i = id + 1; i < ZN; i++){
    
            let tr = document.getElementById(`row${i}`);
            tr.id = `row${i - 1}`;
            let a = tr.querySelector(`#zone_id${i}`);
            a.id = `zone_id${i - 1}`;
            a.value = i;
            let b = tr.querySelector(`#st${i}`);
            b.id = `st${i - 1}`;
            let berr = tr.querySelector(`#sterr${i}`);
            berr.id = `sterr${i - 1}`;
            let c = tr.querySelector(`#ss${i}`);
            c.id = `ss${i - 1}`;
            let cerr = tr.querySelector(`#sserr${i}`);
            cerr.id = `sserr${i - 1}`;
            let save = tr.querySelector(`#save${i}`);
            save.id = `save${i - 1}`;
            save.dataset.id = `${i - 1}`;
            let remove = tr.querySelector(`#remove${i}`);
            remove.id = `remove${i - 1}`;
            remove.dataset.id = `${i - 1}`;
        }
        ZN = ZN - 1;
        
        // RESET ZMAP
        ZMAP = undefined;
            
        // RESET SOLUTION AND RESULTS
        solution_quadrature_reset();
        results_reset();
        if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
        if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
        if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
        document.getElementById("initial_img").style.display = "block";
  
    });
    td2.append(btn2);
    tr.append(td2);
  
    table.append(tr);
}

// FUNCTION TO VERIFY IF ZONE INPUTS ARE VALID
function setup_zone_validator(id){
    valid = true;
  
    // TOTAL CROSS SECTION VALIDATION
    let st_input = document.getElementById(`st${id}`);
    let st_err = document.getElementById(`sterr${id}`);
    let st = st_input.value * 1;
    if (isNaN(st) || st <= 0){
        st_input.classList.remove("is-valid");
        st_input.classList.add("is-invalid");
        st_err.className = "invalid-feedback";
        st_err.innerHTML = "Number > 0";
        valid = false;
    }
    else {
        st_input.classList.remove("is-invalid")
        st_input.classList.add("is-valid");
        st_err.className = "valid-feedback";
        st_err.innerHTML = "";
    }
  
    // SCATTERING CROSS SECTION VALIDATION
    let ss_input = document.getElementById(`ss${id}`);
    let ss_err = document.getElementById(`sserr${id}`);
    let ss = ss_input.value;
    if (ss == ""){ss = "empty";}
    ss = ss * 1;
    if (isNaN(ss) || ss < 0){
        ss_input.classList.remove("is-valid");
        ss_input.classList.add("is-invalid");
        ss_err.className = "invalid-feedback";
        ss_err.innerHTML = "Number >= 0";
        valid = false;
    }
    else if (!isNaN(st) && st > 0){
        if (ss > st){
            ss_input.classList.remove("is-valid");
            ss_input.classList.add("is-invalid");
            ss_err.className = "invalid-feedback";
            ss_err.innerHTML = "Number < total cross section";
            valid = false;
        }
    }
    else {
        ss_input.classList.remove("is-invalid")
        ss_input.classList.add("is-valid");
        ss_err.className = "valid-feedback";
        ss_err.innerHTML = "";
    }
  
    // CLEAR VALIDATORS IF VALID IS TRUE
    if (valid){
        st_input.classList.remove("is-invalid")
        st_input.classList.remove("is-valid");
        ss_input.classList.remove("is-invalid")
        ss_input.classList.remove("is-valid");
    }
  
    return valid;
}

// DONT CLOSE ZONE MODAL UNTIL EVERYTHING IS ALLRIGHT
document.querySelector("#modal-zone").addEventListener("hide.bs.modal", (event) => {
    let msg = document.querySelector("#modal-zone__msg");
    for(let i = 0; i < ZN; i++){
      if (document.querySelector(`#save${i}`).dataset.action == "save"){
        msg.style.display = "block";
        msg.innerHTML = "All rows need to be saved";
        event.preventDefault();
        return;
      }
    }
    msg.style.display = "none";
    msg.innerHTML = "";
});



/////////////////////////////////////////////////////// SETUP SECTION => PHYSICAL REGIONS SETTINGS
/////////////////////////////////////////////////////////////////////////////////////////////////////

// INITIALIZE ANY AVAILABLE REGION
function setup_physical_regions_init() {
    document.querySelector("#modal-physical-regions__x-table").innerHTML = "";
    document.querySelector("#modal-physical-regions__y-table").innerHTML = "";
    for (let i = 0; i < XREG; i++){
        setup_physical_regions_x_load(i, XDOM[i].len, XDOM[i].nc);
    }
    for (let j = 0; j < YREG; j++){
        setup_physical_regions_y_load(j, YDOM[j].len, YDOM[j].nc);
    }
}
  
// ADD A NEW REGION IN X 
function setup_physical_regions_x_add() {
    setup_physical_regions_x_load(XREG, "", "");
    document.querySelector(`#xlen${XREG}`).disabled = false;
    document.querySelector(`#xnc${XREG}`).disabled = false;
    let btn = document.querySelector(`#xr_save${XREG}`);
    btn.dataset.action = "save";
    btn.className = "btn btn-success btn-sm";
    btn.innerHTML = "<i class='bx bx-save' ></i>";
    XREG = XREG + 1;
    let xr = {
        "len": "",
        "nc": ""
    }
    XDOM.push(xr);
}
  
// ADD A NEW REGION IN Y
function setup_physical_regions_y_add() {
    setup_physical_regions_y_load(YREG, "", "");
    document.querySelector(`#ylen${YREG}`).disabled = false;
    document.querySelector(`#ync${YREG}`).disabled = false;
    let btn = document.querySelector(`#yr_save${YREG}`);
    btn.dataset.action = "save";
    btn.className = "btn btn-success btn-sm";
    btn.innerHTML = "<i class='bx bx-save' ></i>";
    YREG = YREG + 1;
    let yr = {
        "len": "",
        "nc": ""
    }
    YDOM.push(yr);
}
  
// LOAD REGIONS IN X
function setup_physical_regions_x_load(xid, xlen, xnc){
  
    let table = document.getElementById("modal-physical-regions__x-table");
    let tr = document.createElement("tr");
    tr.id = `xr_row${xid}`;
  
    // X REGION ID
    let id = document.createElement("td");
    id.className = "input-group-sm";
  
    let id_input = document.createElement("input");
    id_input.id = `xr_id${xid}`;
    id_input.disabled = true;
    id_input.className = "form-control";
    id_input.type = "text";
    id_input.value = `${xid + 1}`;
    id.append(id_input);
    tr.append(id);
  
    // LENGHT OF THE X REGION
    let len = document.createElement("td");
    len.className = "input-group-sm";
    let len_input = document.createElement("input");
    let len_err = document.createElement("div");
    len_input.className = "form-control";
    len_input.id = `xlen${xid}`;
    len_err.id = `xlenerr${xid}`;
    len_input.disabled = true;
    len_input.type = "text";
    len_input.setAttribute("autocomplete", "off");
    len_input.value = xlen;
    len.append(len_input);
    len.append(len_err);
    tr.append(len);
  
    // NUMBER OF CELLS IN THE X REGION
    let nc = document.createElement("td");
    nc.className = "input-group-sm";
    let nc_input = document.createElement("input");
    let nc_err = document.createElement("div");
    nc_input.className = "form-control";
    nc_input.id = `xnc${xid}`;
    nc_err.id = `xncerr${xid}`;
    nc_input.disabled = true;
    nc_input.type = "text";
    nc_input.setAttribute("autocomplete", "off");
    nc_input.value = xnc;
    nc.append(nc_input);
    nc.append(nc_err);
    tr.append(nc);
  
    // EDIT AND SAVE BUTTONS
    let td = document.createElement("td");
    let btn = document.createElement("button");
    btn.id = `xr_save${xid}`;
    btn.className = "btn btn-secondary btn-sm";
    btn.setAttribute("data-id", `${xid}`);
    btn.setAttribute("data-action", "edit");
    btn.innerHTML = "<i class='bx bx-edit' ></i>";
    btn.onclick = function(){
        if (setup_physical_regions_x_validator(this.dataset.id) === true){
            let id = parseInt(this.dataset.id);
            let len = document.getElementById(`xlen${id}`);
            let nc = document.getElementById(`xnc${id}`);
            let action = this.dataset.action;
            if (action == "edit"){
                len.disabled = false;
                nc.disabled = false;
                this.dataset.action = "save";
                this.className = "btn btn-success btn-sm";
                this.innerHTML = "<i class='bx bx-save' ></i>";
            }
            else {
                XDOM[id].len = parseFloat(len.value);
                XDOM[id].nc = parseInt(nc.value);
                len.disabled = true;
                nc.disabled = true;
                this.dataset.action = "edit";
                this.className = "btn btn-secondary btn-sm";
                this.innerHTML = "<i class='bx bx-edit'></i>";
            }
            
            // RESET ZMAP AND QMAP
            ZMAP = undefined; QMAP = undefined;
            
            // RESET SOLUTION AND RESULTS
            solution_quadrature_reset();
            results_reset();
            if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
            if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
            if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
            document.getElementById("initial_img").style.display = "block";
        }
    };
    td.append(btn);
    tr.append(td);
  
    // REMOVE BUTTON
    let td2 = document.createElement("td");
    let btn2 = document.createElement("button");
    btn2.id = `xr_remove${xid}`;
    btn2.className = "btn btn-danger btn-sm";
    btn2.setAttribute("data-id", `${xid}`);
    btn2.innerHTML = "<i class='bx bxs-trash'></i>";
    btn2.addEventListener("click", function(){
  
        let id = parseInt(this.dataset.id);
        document.getElementById(`xr_row${id}`).remove();
        XDOM.splice(id, 1);
        for(let i = id + 1; i < XREG; i++){
    
            let tr = document.getElementById(`xr_row${i}`);
            tr.id = `xr_row${i - 1}`;
            let a = tr.querySelector(`#xr_id${i}`);
            a.id = `xr_id${i - 1}`;
            a.value = i;
            let b = tr.querySelector(`#xlen${i}`);
            b.id = `xlen${i - 1}`;
            let berr = tr.querySelector(`#xlenerr${i}`);
            berr.id = `xlenerr${i - 1}`;
            let c = tr.querySelector(`#xnc${i}`);
            c.id = `xnc${i - 1}`;
            let cerr = tr.querySelector(`#xncerr${i}`);
            cerr.id = `xncerr${i - 1}`;
            let save = tr.querySelector(`#xr_save${i}`);
            save.id = `xr_save${i - 1}`;
            save.dataset.id = `${i - 1}`;
            let remove = tr.querySelector(`#xr_remove${i}`);
            remove.id = `xr_remove${i - 1}`;
            remove.dataset.id = `${i - 1}`;
        }
        XREG = XREG - 1;
        
        // RESET ZMAP AND QMAP
        ZMAP = undefined; QMAP = undefined;
            
        // RESET SOLUTION AND RESULTS
        solution_quadrature_reset();
        results_reset();
        if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
        if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
        if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
        document.getElementById("initial_img").style.display = "block";
  
    });
    td2.append(btn2);
    tr.append(td2);
  
    table.append(tr);
}
  
// LOAD REGIONS IN Y
function setup_physical_regions_y_load(yid, ylen, ync){
  
    let table = document.getElementById("modal-physical-regions__y-table");
    let tr = document.createElement("tr");
    tr.id = `yr_row${yid}`;
  
    // Y REGION ID
    let id = document.createElement("td");
    id.className = "input-group-sm";
  
    let id_input = document.createElement("input");
    id_input.id = `yr_id${yid}`;
    id_input.disabled = true;
    id_input.className = "form-control";
    id_input.type = "text";
    id_input.value = `${yid + 1}`;
    id.append(id_input);
    tr.append(id);
  
    // LENGHT OF THE Y REGION
    let len = document.createElement("td");
    len.className = "input-group-sm";
    let len_input = document.createElement("input");
    let len_err = document.createElement("div");
    len_input.className = "form-control";
    len_input.id = `ylen${yid}`;
    len_err.id = `ylenerr${yid}`;
    len_input.disabled = true;
    len_input.type = "text";
    len_input.setAttribute("autocomplete", "off");
    len_input.value = ylen;
    len.append(len_input);
    len.append(len_err);
    tr.append(len);
  
    // NUMBER OF CELLS IN THE Y REGION
    let nc = document.createElement("td");
    nc.className = "input-group-sm";
    let nc_input = document.createElement("input");
    let nc_err = document.createElement("div");
    nc_input.className = "form-control";
    nc_input.id = `ync${yid}`;
    nc_err.id = `yncerr${yid}`;
    nc_input.disabled = true;
    nc_input.type = "text";
    nc_input.setAttribute("autocomplete", "off");
    nc_input.value = ync;
    nc.append(nc_input);
    nc.append(nc_err);
    tr.append(nc);
  
    // EDIT AND SAVE BUTTONS
    let td = document.createElement("td");
    let btn = document.createElement("button");
    btn.id = `yr_save${yid}`;
    btn.className = "btn btn-secondary btn-sm";
    btn.setAttribute("data-id", `${yid}`);
    btn.setAttribute("data-action", "edit");
    btn.innerHTML = "<i class='bx bx-edit' ></i>";
    btn.onclick = function(){
        if (setup_physical_regions_y_validator(this.dataset.id) === true){
            let id = parseInt(this.dataset.id);
            let len = document.getElementById(`ylen${id}`);
            let nc = document.getElementById(`ync${id}`);
            let action = this.dataset.action;
            if (action == "edit"){
                len.disabled = false;
                nc.disabled = false;
                this.dataset.action = "save";
                this.className = "btn btn-success btn-sm";
                this.innerHTML = "<i class='bx bx-save' ></i>";
            }
            else {
                YDOM[id].len = parseFloat(len.value);
                YDOM[id].nc = parseInt(nc.value);
                len.disabled = true;
                nc.disabled = true;
                this.dataset.action = "edit";
                this.className = "btn btn-secondary btn-sm";
                this.innerHTML = "<i class='bx bx-edit'></i>";
            }
            
            // RESET ZMAP AND QMAP
            ZMAP = undefined; QMAP = undefined;
            
            // RESET SOLUTION AND RESULTS
            solution_quadrature_reset();
            results_reset();
            if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
            if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
            if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
            document.getElementById("initial_img").style.display = "block";
        }
    };
    td.append(btn);
    tr.append(td);
  
    // REMOVE BUTTON
    let td2 = document.createElement("td");
    let btn2 = document.createElement("button");
    btn2.id = `yr_remove${yid}`;
    btn2.className = "btn btn-danger btn-sm";
    btn2.setAttribute("data-id", `${yid}`);
    btn2.innerHTML = "<i class='bx bxs-trash'></i>";
    btn2.addEventListener("click", function(){
  
        let id = parseInt(this.dataset.id);
        document.getElementById(`yr_row${id}`).remove();
        YDOM.splice(id, 1);
        for(let i = id + 1; i < YREG; i++){
    
            let tr = document.getElementById(`yr_row${i}`);
            tr.id = `yr_row${i - 1}`;
            let a = tr.querySelector(`#yr_id${i}`);
            a.id = `yr_id${i - 1}`;
            a.value = i;
            let b = tr.querySelector(`#ylen${i}`);
            b.id = `ylen${i - 1}`;
            let berr = tr.querySelector(`#ylenerr${i}`);
            berr.id = `ylenerr${i - 1}`;
            let c = tr.querySelector(`#ync${i}`);
            c.id = `ync${i - 1}`;
            let cerr = tr.querySelector(`#yncerr${i}`);
            cerr.id = `yncerr${i - 1}`;
            let save = tr.querySelector(`#yr_save${i}`);
            save.id = `yr_save${i - 1}`;
            save.dataset.id = `${i - 1}`;
            let remove = tr.querySelector(`#yr_remove${i}`);
            remove.id = `yr_remove${i - 1}`;
            remove.dataset.id = `${i - 1}`;
        }
        YREG = YREG - 1;
        
        // RESET ZMAP AND QMAP
        ZMAP = undefined; QMAP = undefined;
            
        // RESET SOLUTION AND RESULTS
        solution_quadrature_reset();
        results_reset();
        if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
        if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
        if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
        document.getElementById("initial_img").style.display = "block";
    
    });
    td2.append(btn2);
    tr.append(td2);
  
    table.append(tr);
}
  
// VALIDATE X REGION FUNCTION
function setup_physical_regions_x_validator(id){
    valid = true;
  
    // LENGHT VALIDATION
    let len_input = document.getElementById(`xlen${id}`);
    let len_err = document.getElementById(`xlenerr${id}`);
    let len = len_input.value * 1;
    if (isNaN(len) || len <= 0){
        len_input.classList.remove("is-valid");
        len_input.classList.add("is-invalid");
        len_err.className = "invalid-feedback";
        len_err.innerHTML = "Number > 0";
        valid = false;
    }
    else {
        len_input.classList.remove("is-invalid")
        len_input.classList.add("is-valid");
        len_err.className = "valid-feedback";
        len_err.innerHTML = "";
    }
  
    // NUMBER OF CELLS VALIDATION
    let nc_input = document.getElementById(`xnc${id}`);
    let nc_err = document.getElementById(`xncerr${id}`);
    let nc = nc_input.value * 1;
    if (!Number.isInteger(nc) || nc <= 0){
        nc_input.classList.remove("is-valid");
        nc_input.classList.add("is-invalid");
        nc_err.className = "invalid-feedback";
        nc_err.innerHTML = "Integer > 0";
        valid = false;
    }
    else {
        nc_input.classList.remove("is-invalid")
        nc_input.classList.add("is-valid");
        nc_err.className = "valid-feedback";
        nc_err.innerHTML = "";
    }
  
    // CLEAR VALIDATORS IS VALID IS TRUE
    if (valid){
        len_input.classList.remove("is-invalid")
        len_input.classList.remove("is-valid");
        nc_input.classList.remove("is-invalid")
        nc_input.classList.remove("is-valid");
    }
  
    return valid;
}
  
// VALIDATE Y REGION FUNCTION
function setup_physical_regions_y_validator(id){
    valid = true;
  
    // LENGHT VALIDATION
    let len_input = document.getElementById(`ylen${id}`);
    let len_err = document.getElementById(`ylenerr${id}`);
    let len = len_input.value * 1;
    if (isNaN(len) || len <= 0){
        len_input.classList.remove("is-valid");
        len_input.classList.add("is-invalid");
        len_err.className = "invalid-feedback";
        len_err.innerHTML = "Number > 0";
        valid = false;
    }
    else {
        len_input.classList.remove("is-invalid")
        len_input.classList.add("is-valid");
        len_err.className = "valid-feedback";
        len_err.innerHTML = "";
    }
  
    // NUMBER OF CELLS VALIDATION
    let nc_input = document.getElementById(`ync${id}`);
    let nc_err = document.getElementById(`yncerr${id}`);
    let nc = nc_input.value * 1;
    if (!Number.isInteger(nc) || nc <= 0){
        nc_input.classList.remove("is-valid");
        nc_input.classList.add("is-invalid");
        nc_err.className = "invalid-feedback";
        nc_err.innerHTML = "Integer > 0";
        valid = false;
    }
    else {
        nc_input.classList.remove("is-invalid")
        nc_input.classList.add("is-valid");
        nc_err.className = "valid-feedback";
        nc_err.innerHTML = "";
    }
  
    // CLEAR VALIDATORS IS VALID IS TRUE
    if (valid){
        len_input.classList.remove("is-invalid")
        len_input.classList.remove("is-valid");
        nc_input.classList.remove("is-invalid")
        nc_input.classList.remove("is-valid");
    }
  
    return valid;
}
  
// DONT CLOSE REGION MODAL UNTIL EVERYTHING IS ALLRIGHT
document.querySelector("#modal-physical-regions").addEventListener("hide.bs.modal", (event) => {
    let msg = document.querySelector("#modal-physical-regions__msg");
    for(let i = 0; i < XREG; i++){
        if (document.querySelector(`#xr_save${i}`).dataset.action == "save"){
            msg.style.display = "block";
            msg.innerHTML = "All rows need to be saved";
            event.preventDefault();
            return;
        }
    }
    for(let i = 0; i < YREG; i++){
        if (document.querySelector(`#yr_save${i}`).dataset.action == "save"){
            msg.style.display = "block";
            msg.innerHTML = "All rows need to be saved.";
            event.preventDefault();
            return;
        }
    }
    msg.style.display = "none";
    msg.innerHTML = "";
});



////////////////////////////////////////////////////////////// SETUP SECTION => ZONE MAPPING SETTINGS
/////////////////////////////////////////////////////////////////////////////////////////////////////

// CREATE ZMAP FUNCTION
function setup_zone_map_create(){
    ZMAP = [];
    for(let j = 0; j < YREG; j++){
        ZMAP[j] = [];
        for(let i = 0; i < XREG; i++){
            ZMAP[j][i] = ZON.length - 1;
        }
    }
}
  
// ZONE MAPPING FUNCTION
function setup_zone_map_init(){
  
    // VALIDATION
    let msg = document.querySelector("#modal-zone-map__msg");
    if (XREG == 0) {
        msg.style.display = "block";
        msg.innerHTML = "There are no regions in x specified yet.";
        return;
    }
    else if (YREG == 0) {
        msg.style.display = "block";
        msg.innerHTML = "There are no regions in y specified yet.";
        return;
    }
    else if (ZN == 0) {
        msg.style.display = "block";
        msg.innerHTML = "There are no zone specified yet.";
        return;
    }
    else {
        msg.style.display = "none";
        msg.innerHTML = "";
    }
  
    // INITIALIZE ZMAP IF THERE IS ANY ERROR
    if (ZMAP == undefined) setup_zone_map_create();
    if (ZMAP.length != YREG) setup_zone_map_create();
    if (ZMAP[0].length != XREG) setup_zone_map_create();
  
    // CREATE TABLE
    let table = document.getElementById("modal-zone-map__table");
    table.innerHTML = "";
    for(let j = YREG; j >= 0; j--){
        let tr = document.createElement("tr");
        tr.className = "text-nowrap";
        for(let i = 0; i <= XREG; i++){
            let td = document.createElement("td");
            
            if(i == 0 && j == 0){
                td.innerHTML = "RY \\ RX";
                td.className = "px-2";
            }
            else if (j == 0){
                td.innerHTML = `${i}`;
                td.className = "px-2";
            }
            else if (i == 0){
                td.innerHTML = `${j}`;
                td.className = "px-2";
            }
            else {
                if (ZMAP[j-1][i-1] == undefined){
                    ZMAP[j-1][i-1] = ZON.length - 1;
                }
                let select = document.createElement("select");
                select.className = "form-select form-select-sm border-0";
                select.setAttribute("data-xr", `${i-1}`);
                select.setAttribute("data-yr", `${j-1}`);
                for(let z = 0; z < ZN; z++){
                    let option = document.createElement("option");
                    option.text = `Zone ${z + 1}`;
                    if (z == ZMAP[j-1][i-1]){
                        option.selected = true;
                    }
                    select.append(option);
                }
                select.onchange = function(){
                    let i = this.dataset.xr;
                    let j = this.dataset.yr;
                    let z = this.selectedIndex;
                    ZMAP[j][i] = z;
                    
                    // RESET SOLUTION AND RESULTS
                    solution_quadrature_reset();
                    results_reset();
                    if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
                    if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
                    if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
                    document.getElementById("initial_img").style.display = "block";

                }
                td.append(select);
            }
            tr.append(td);
        }
        table.append(tr);
    }
}



//////////////////////////////////////////////////////////// SETUP SECTION => SOURCE MAPPING SETTINGS
/////////////////////////////////////////////////////////////////////////////////////////////////////

// CREATE QMAP FUNCTION
function setup_source_map_create(){
    QMAP = [];
    for(let j = 0; j < YREG; j++){
        QMAP[j] = [];
        for(let i = 0; i < XREG; i++){
            QMAP[j][i] = 0;
        }
    }
}
  
// SOURCE MAPPING FUNCTION
function setup_source_map_init(){
  
    let msg = document.querySelector("#modal-source-map__msg");
    if (XREG == 0) {
        msg.style.display = "block";
        msg.innerHTML = "There are no regions in x specified yet";
        return;
    }
    else if (YREG == 0) {
        msg.style.display = "block";
        msg.innerHTML = "There are no regions in y specified yet";
        return;
    }
    else {
        msg.style.display = "none";
        msg.innerHTML = "";
    }
  
    // INICIALIZE QMAP IF ANY PROBLEM
    if (QMAP == undefined) setup_source_map_create();
    if (QMAP.length != YREG) setup_source_map_create();
    if (QMAP[0].length != XREG) setup_source_map_create();
  
    // CREATE TABLE
    let table = document.getElementById("modal-source-map__table");
    table.innerHTML = "";
    for(let j = YREG; j >= 0; j--){
        let tr = document.createElement("tr");
        for(let i = 0; i <= XREG; i++){
            let td = document.createElement("td");
            td.className = "input-group-sm";
            if(i == 0 && j == 0){
                td.innerHTML = "RY \\ RX";
            }
            else if (j == 0){
                td.innerHTML = `${i}`;
            }
            else if (i == 0){
                td.innerHTML = `${j}`;
            }
            else {
                if (QMAP[j-1][i-1] == undefined){
                    QMAP[j-1][i-1] = 0;
                }
                let input = document.createElement("input");
                let qerr = document.createElement("div");
                input.className = "form-control";
                qerr.className = "invalid-feedback";
                qerr.id = `qerr${j-1}${i-1}`;
                input.type = "text";
                input.setAttribute("autocomplete", "off");
                input.setAttribute("data-xr", `${i-1}`);
                input.setAttribute("data-yr", `${j-1}`);
                input.value = `${QMAP[j-1][i-1]}`;
                input.onchange = function(){
                    let i = this.dataset.xr;
                    let j = this.dataset.yr;
                    if (setup_source_map_validator(this) == true) QMAP[j][i] = parseFloat(this.value);
                    else QMAP[j][i] = 0;
                    
                    // RESET SOLUTION AND RESULTS
                    solution_quadrature_reset();
                    results_reset();
                    if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
                    if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
                    if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
                    document.getElementById("initial_img").style.display = "block";

                }
                td.append(input);
                td.append(qerr);
            }
            tr.append(td);
        }
        table.append(tr);
    }
}
  
// QMAP VALIDATOR FUNCTION
function setup_source_map_validator(input){
    let valid = true;
    let i = input.dataset.xr;
    let j = input.dataset.yr;
    let err = document.querySelector(`#qerr${j}${i}`);
    let q = input.value;
    if (q == "") q = "empty";
    q = q * 1;
    if (isNaN(q) || q < 0){
        input.classList.add("is-invalid");
        err.innerHTML = "Number >= 0";
        valid = false;
    }
    else {
        input.classList.remove("is-invalid")
        err.innerHTML = "";
    }
    return valid;
}



//////////////////////////////////////////////////////// SETUP SECTION => BOUNDARY CONDITION SETTINGS
/////////////////////////////////////////////////////////////////////////////////////////////////////

// BOUNDARY CONDITIONS INITIALIZATION
function setup_bc_init() {

    let left = document.querySelector("#bc_left_select");
    if (BC["left"] == 0) left.options[0].selected = true;
    else if (BC["left"] == -1) left.options[1].selected = true;
    else if (BC["left"] > 0) left.options[2].selected = true;
    else left.options[0].selected = true;
  
    let right = document.querySelector("#bc_right_select");
    if (BC["right"] == 0) right.options[0].selected = true;
    else if (BC["right"] == -1) right.options[1].selected = true;
    else if (BC["right"] > 0) right.options[2].selected = true;
    else right.options[0].selected = true;
  
    let top = document.querySelector("#bc_top_select");
    if (BC["top"] == 0) top.options[0].selected = true;
    else if (BC["top"] == -1) top.options[1].selected = true;
    else if (BC["top"] > 0) top.options[2].selected = true;
    else top.options[0].selected = true;
  
    let bottom = document.querySelector("#bc_bottom_select");
    if (BC["bottom"] == 0) bottom.options[0].selected = true;
    else if (BC["bottom"] == -1) bottom.options[1].selected = true;
    else if (BC["bottom"] > 0) bottom.options[2].selected = true;
    else bottom.options[0].selected = true;
}
  
// LEFT
document.querySelector("#bc_left_select").onchange = function (){
    let left = document.querySelector("#bc_left_input");
    let err = document.querySelector("#bc_left_err");
    left.classList.remove("is-invalid");
    err.innerHTML = "";
    if (this.selectedIndex == 0){
        left.style.display = "none";
        BC["left"] = 0;
    }
    else if (this.selectedIndex == 1){
        left.style.display = "none";
        BC["left"] = -1;
    }
    else {
        left.style.display = "inline-block";
        BC["left"] = 0;
        left.value = `${BC["left"]}`;
    }
    
    // RESET SOLUTION AND RESULTS
    solution_quadrature_reset();
    results_reset();
    if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
    if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
    if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
    document.getElementById("initial_img").style.display = "block";
}
document.querySelector("#bc_left_input").onchange = function (){
    if (setup_bc_validator("left")) BC["left"] = parseFloat(this.value);
}
  
// RIGHT
document.querySelector("#bc_right_select").onchange = function (){
    let right = document.querySelector("#bc_right_input");
    let err = document.querySelector("#bc_right_err");
    right.classList.remove("is-invalid");
    err.innerHTML = "";
    if (this.selectedIndex == 0){
        right.style.display = "none";
        BC["right"] = 0;
    }
    else if (this.selectedIndex == 1){
        right.style.display = "none";
        BC["right"] = -1;
    }
    else {
        right.style.display = "inline-block";
        BC["right"] = 0;
        right.value = `${BC["right"]}`;
    }
    
    // RESET SOLUTION AND RESULTS
    solution_quadrature_reset();
    results_reset();
    if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
    if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
    if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
    document.getElementById("initial_img").style.display = "block";
}
document.querySelector("#bc_right_input").onchange = function (){
    if (setup_bc_validator("right")) BC["right"] = parseFloat(this.value);
}
  
// TOP
document.querySelector("#bc_top_select").onchange = function (){
    let top = document.querySelector("#bc_top_input");
    let err = document.querySelector("#bc_top_err");
    top.classList.remove("is-invalid");
    err.innerHTML = "";
    if (this.selectedIndex == 0){
        top.style.display = "none";
        BC["top"] = 0;
    }
    else if (this.selectedIndex == 1){
        top.style.display = "none";
        BC["top"] = -1;
    }
    else {
        top.style.display = "inline-block";
        BC["top"] = 0;
        top.value = `${BC["top"]}`;
    }
    
    // RESET SOLUTION AND RESULTS
    solution_quadrature_reset();
    results_reset();
    if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
    if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
    if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
    document.getElementById("initial_img").style.display = "block";
}
document.querySelector("#bc_top_input").onchange = function (){
    if (setup_bc_validator("top")) BC["top"] = parseFloat(this.value);
}
  
// BOTTOM
document.querySelector("#bc_bottom_select").onchange = function (){
    let bottom = document.querySelector("#bc_bottom_input");
    let err = document.querySelector("#bc_bottom_err");
    bottom.classList.remove("is-invalid");
    err.innerHTML = "";
    if (this.selectedIndex == 0){
        bottom.style.display = "none";
        BC["bottom"] = 0;
    }
    else if (this.selectedIndex == 1){
        bottom.style.display = "none";
        BC["bottom"] = -1;
    }
    else {
        bottom.style.display = "inline-block";
        BC["bottom"] = 0;
        bottom.value = `${BC["bottom"]}`;
    }
    
    // RESET SOLUTION AND RESULTS
    solution_quadrature_reset();
    results_reset();
    if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
    if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
    if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
    document.getElementById("initial_img").style.display = "block";
}
document.querySelector("#bc_bottom_input").onchange = function (){
    if (setup_bc_validator("bottom")) BC["bottom"] = parseFloat(this.value);
}
  
// BC VALIDATOR FUNCTION
function setup_bc_validator(text){
    let valid = true;
    let input = document.querySelector(`#bc_${text}_input`);
    let err = document.querySelector(`#bc_${text}_err`);
    let bc = input.value;
    input.classList.remove("is-invalid")
    err.innerHTML = "";
    if (bc == "") bc = "empty";
    bc = bc * 1;
    if (isNaN(bc) || bc < 0){
        input.classList.add("is-invalid");
        err.innerHTML = "Number >= 0";
        valid = false;
    }
    return valid;
}



/////////////////////////////////////////////////////////////////// SETUP SECTION => DISPLAY GEOMETRY
/////////////////////////////////////////////////////////////////////////////////////////////////////

// DISPLAY GEOMETRY FUNCTION
function setup_display_geometry(){

    // WARNINGS
    if (ZN == 0 || XREG == 0 || YREG == 0 || ZMAP == undefined){
        console_process_error("geometry");
        return;
    }

    document.getElementById("initial_img").style.display = "none";
    if (document.getElementById("canvas_geom")) document.getElementById("canvas_geom").remove();
    if (document.getElementById("canvas_quad")) document.getElementById("canvas_quad").remove();
    if (document.getElementById("results__picture")) document.getElementById("results__picture").remove();
  
    // CREATE CANVAS
    let monitor = document.getElementById("main-content__layout--monitor");
    let canvas = document.createElement("canvas");
    canvas.id = "canvas_geom";
    canvas.width = monitor.offsetWidth - 50;
    canvas.height = monitor.offsetHeight - 50;
  
    // CANVAS SETTINGS
    let XLEN = 0;
    for(let i = 0; i < XREG; i++){
      XLEN = XLEN + XDOM[i].len;
    }
    let YLEN = 0;
    for(let i = 0; i < YREG; i++){
      YLEN = YLEN + YDOM[i].len;
    }
    let SCALE = Math.min((canvas.height-50)/YLEN, (canvas.width-50)/XLEN);
    let X_MARGIN = (canvas.width - XLEN * SCALE) / 2;
    let Y_MARGIN = (canvas.height - YLEN * SCALE) / 2;

    let myFont = 14, mySpan = 15;
    if (canvas.width < 500) {
        myFont = 10;
        mySpan = 10;
    }
  
    // CANVAS CONTENT
    if (canvas.getContext) {
        let ctx = canvas.getContext("2d");
        ctx.clearRect(0, 0, canvas.width, canvas.height);
    
        // DOMAIN
        ctx.shadowOffsetX = 4;
        ctx.shadowOffsetY = 4;
        ctx.shadowBlur = 4;
        ctx.shadowColor = 'rgba(0, 0, 0, 0.6)';
        ctx.fillStyle = "#343a40";
        ctx.fillRect(X_MARGIN, Y_MARGIN, XLEN * SCALE, YLEN * SCALE);
        ctx.shadowOffsetX = 0;
        ctx.shadowOffsetY = 0;
        ctx.shadowBlur = 0;
        ctx.strokeStyle = "#212529";
        ctx.lineWidth = 2;
        ctx.strokeRect(X_MARGIN, Y_MARGIN, XLEN * SCALE, YLEN * SCALE);

        // WRITE DIMENSIONS INIT
        ctx.font = `${myFont} Times New Roman`;
        ctx.fillStyle = "Black";
        ctx.fillText(`0`, X_MARGIN, Y_MARGIN - mySpan);
        ctx.fillText(`0`, X_MARGIN + XLEN * SCALE + mySpan, Y_MARGIN + YLEN * SCALE);
    
        // REGIONS
        let x0, y0, xsum = 0, ysum = 0;
        x0 = X_MARGIN;
        for (let i = 0; i < XREG; i++) {
            y0 = canvas.height - Y_MARGIN;
            for (let j = 0; j < YREG; j++) {
                y0 = y0 - YDOM[j].len * SCALE;
                ctx.fillStyle = `rgba(250,250,250,${1 - ZON[ZMAP[j][i]].ss/ZON[ZMAP[j][i]].st})`;
                if (document.querySelector("#show_source").checked){
                    if (QMAP[j][i] != 0){
                        ctx.fillStyle = "#dc3545";
                    }
                }
                ctx.fillRect(x0, y0, XDOM[i].len * SCALE, YDOM[j].len * SCALE);
                ctx.strokeRect(x0, y0, XDOM[i].len * SCALE, YDOM[j].len * SCALE);

                // WRITE DIMENSIONS
                if (j == (YREG - 1)){
                    xsum = xsum + XDOM[i].len;
                    ctx.font = `${myFont} Times New Roman`;
                    ctx.fillStyle = "Black";
                    ctx.fillText(`${xsum}`, x0 + XDOM[i].len * SCALE, y0 - mySpan);
                }
                if (i == (XREG - 1)) {
                    ysum = ysum + YDOM[j].len;
                    ctx.font = `${myFont} Times New Roman`;
                    ctx.fillStyle = "Black";
                    ctx.fillText(`${ysum}`, x0 + XDOM[i].len * SCALE + mySpan, y0);
                }
            }
            x0 = x0 + XDOM[i].len * SCALE;
        }

        // SHOW MESH
        if (document.querySelector("#show_mesh").checked){
            ctx.lineWidth = 1;
            x0 = X_MARGIN; y0 = Y_MARGIN;
            for (let i = 0; i < XREG; i++) {
            let hi = XDOM[i].len / XDOM[i].nc;
            for (let c = 0; c < XDOM[i].nc; c++){
                x0 = x0 + hi * SCALE;
                ctx.beginPath();
                ctx.moveTo(x0, y0);
                ctx.lineTo(x0, y0 + YLEN * SCALE);
                ctx.stroke();
            }
            }
            x0 = X_MARGIN; y0 = Y_MARGIN + YLEN * SCALE;
            for (let j = 0; j < YREG; j++) {
            let hj = YDOM[j].len / YDOM[j].nc;
            for (let c = 0; c < YDOM[j].nc; c++){
                y0 = y0 - hj * SCALE;
                ctx.beginPath();
                ctx.moveTo(x0, y0);
                ctx.lineTo(x0 + XLEN * SCALE, y0);
                ctx.stroke();
            }
            }
        }
    
        // SHOW SOURCE
        if (document.querySelector("#show_source").checked){
            ctx.lineWidth = 2;
            ctx.strokeStyle = "red";
            if (BC["left"] > 0){
            ctx.beginPath();
            ctx.moveTo(X_MARGIN, Y_MARGIN);
            ctx.lineTo(X_MARGIN, Y_MARGIN + YLEN * SCALE);
            ctx.stroke();
            }
            if (BC["right"] > 0){
            ctx.beginPath();
            ctx.moveTo(X_MARGIN + XLEN * SCALE, Y_MARGIN);
            ctx.lineTo(X_MARGIN + XLEN * SCALE, Y_MARGIN + YLEN * SCALE);
            ctx.stroke();
            }
            if (BC["top"] > 0){
            ctx.beginPath();
            ctx.moveTo(X_MARGIN, Y_MARGIN);
            ctx.lineTo(X_MARGIN + XLEN * SCALE, Y_MARGIN);
            ctx.stroke();
            }
            if (BC["bottom"] > 0){
            ctx.beginPath();
            ctx.moveTo(X_MARGIN, Y_MARGIN + YLEN * SCALE);
            ctx.lineTo(X_MARGIN + XLEN * SCALE, Y_MARGIN + YLEN * SCALE);
            ctx.stroke();
            }
        }

        // COORDINATE SYSTEM
        ctx.fillStyle = "#343a40";
        ctx.strokeStyle = "#343a40";
        let xlen, ylen;
        ctx.beginPath();
        x0 = X_MARGIN/2 - 10; xlen = 30;
        y0 = Y_MARGIN + YLEN * SCALE + Y_MARGIN/2;
        ctx.moveTo(x0, y0);
        ctx.lineTo(x0 + xlen, y0);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(x0 + xlen + 5, y0);
        ctx.lineTo(x0 + xlen, y0 + 2);
        ctx.lineTo(x0 + xlen, y0 -2);
        ctx.fill();
        ctx.font = `${myFont} Times New Roman`;
        ctx.fillStyle = "Black";
        ctx.fillText('X (cm)', x0 + xlen + 10, y0 + 3);

        ctx.beginPath();
        x0 = X_MARGIN/2; ylen = 30;
        y0 = Y_MARGIN + YLEN * SCALE + Y_MARGIN/2 + 10;
        ctx.moveTo(x0, y0);
        ctx.lineTo(x0, y0 - ylen);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(x0, y0 - ylen - 5);
        ctx.lineTo(x0 + 2, y0 - ylen);
        ctx.lineTo(x0 - 2, y0 - ylen);
        ctx.fill();
        ctx.font = `${myFont} Times New Roman`;
        ctx.fillStyle = "Black";
        ctx.fillText('Y', x0 - 3, y0 - ylen - 10);
    }

    monitor.append(canvas);  
}