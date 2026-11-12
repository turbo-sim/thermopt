(function(){
  "use strict";
  const payload=JSON.parse(document.getElementById("report-payload").textContent);
  const app=document.getElementById("report-app");
  const urls=new Map();
  function assetURL(key){
    if(!urls.has(key)){
      const file=payload.assets[key];
      if(!file)throw new Error("Missing embedded asset: "+key);
      let content=file.text;
      let mime=file.mime+";charset=utf-8";
      if(file.encoding==="base64"){
        content=Uint8Array.from(atob(file.data),character=>character.charCodeAt(0));
        mime=file.mime;
      }
      urls.set(key,URL.createObjectURL(new Blob([content],{type:mime})));
    }
    return urls.get(key);
  }
  function saveFile(key,name){
    const link=document.createElement("a");
    link.href=assetURL(key);link.download=name;
    document.body.append(link);link.click();link.remove();
  }
  function render(){
    const route=new URLSearchParams(location.hash.slice(1));
    const key=Object.hasOwn(payload.pages,route.get("page"))?route.get("page"):"index";
    const page=payload.pages[key];
    const template=document.createElement("template");
    template.innerHTML=page.body;
    template.content.querySelectorAll("[data-embedded-image]").forEach(img=>{
      img.src=assetURL(img.dataset.embeddedImage);
    });
    template.content.querySelectorAll("[data-embedded-file]").forEach(link=>{
      link.addEventListener("click",event=>{
        event.preventDefault();saveFile(link.dataset.embeddedFile,link.download);
      });
    });
    app.replaceChildren(template.content);
    document.body.classList.remove("menu-open");
    document.title=page.title+" · EMPOWER";
    const anchor=route.get("anchor")||"";
    window.initializeReport(anchor);
    document.querySelectorAll("[data-focus-main]").forEach(link=>link.addEventListener("click",event=>{
      event.preventDefault();const main=document.getElementById("main");main.tabIndex=-1;main.focus();
    }));
    document.querySelectorAll("[data-print-page]").forEach(button=>button.addEventListener("click",async()=>{
      const images=[...app.querySelectorAll("img[src]")];
      images.forEach(img=>img.loading="eager");
      await Promise.allSettled(images.map(img=>img.decode?img.decode():Promise.resolve()));
      window.print();
    }));
    if(anchor){document.getElementById(anchor)?.scrollIntoView();}
    else{window.scrollTo(0,0);}
  }
  window.addEventListener("hashchange",render);
  render();
})();
