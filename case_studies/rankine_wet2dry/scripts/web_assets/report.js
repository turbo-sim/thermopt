"use strict";
window.initializeReport=function(anchor=""){
document.querySelectorAll("[data-tabs]").forEach(group=>{
  const buttons=[...group.querySelectorAll("[role=tab]")];
  const activate=button=>{
    buttons.forEach(item=>{
      const selected=item===button;
      item.setAttribute("aria-selected",String(selected));
      item.tabIndex=selected?0:-1;
      document.getElementById(item.getAttribute("aria-controls")).hidden=!selected;
    });
  };
  buttons.forEach((button,index)=>{
    button.addEventListener("click",()=>activate(button));
    button.addEventListener("keydown",event=>{
      let next=index;
      if(event.key==="ArrowRight")next=(index+1)%buttons.length;
      else if(event.key==="ArrowLeft")next=(index+buttons.length-1)%buttons.length;
      else if(event.key==="Home")next=0;
      else if(event.key==="End")next=buttons.length-1;
      else return;
      event.preventDefault();activate(buttons[next]);buttons[next].focus();
    });
  });
  const linked=buttons.find(b=>anchor===b.getAttribute("aria-controls")||"#"+b.getAttribute("aria-controls")===location.hash);
  if(linked)activate(linked);
});
const menu=document.querySelector(".menu-toggle");
const closeMenu=()=>{document.body.classList.remove("menu-open");menu?.setAttribute("aria-expanded","false");};
menu?.addEventListener("click",()=>{const open=document.body.classList.toggle("menu-open");menu.setAttribute("aria-expanded",String(open));});
document.querySelector(".mobile-scrim")?.addEventListener("click",closeMenu);
const dialog=document.getElementById("figure-dialog");
document.querySelectorAll("[data-zoom]").forEach(button=>button.addEventListener("click",()=>{
  const img=button.closest("figure").querySelector("img");
  dialog.querySelector("img").src=img.src;
  dialog.querySelector("img").alt=img.alt;
  dialog.querySelector("#zoom-title").textContent=img.alt;
  dialog.querySelector(".original-link").href=img.src;
  dialog.showModal();
}));
dialog?.querySelector(".zoom-close").addEventListener("click",()=>dialog.close());
dialog?.addEventListener("click",event=>{if(event.target===dialog)dialog.close();});
const controls=[...document.querySelectorAll("[data-filter]")];
const rows=[...document.querySelectorAll("#case-table tbody tr")];
function filterCases(){
  let count=0;
  rows.forEach(row=>{
    const visible=controls.every(input=>{
      const value=input.value.toLowerCase().trim();
      if(!value)return true;
      const field=input.dataset.filter;
      return field==="search"?((row.dataset.search||"")+" "+row.textContent).toLowerCase().includes(value):row.dataset[field]===value;
    });
    row.hidden=!visible;if(visible)count++;
  });
  const label=document.getElementById("case-count");
  if(label)label.textContent=`${count} of ${rows.length} optimized cases shown`;
  const empty=document.getElementById("no-cases");if(empty)empty.hidden=count!==0;
}
controls.forEach(input=>{input.addEventListener("input",filterCases);input.addEventListener("change",filterCases);});
if(rows.length)filterCases();
};
if(!window.REPORT_STANDALONE)window.initializeReport();
