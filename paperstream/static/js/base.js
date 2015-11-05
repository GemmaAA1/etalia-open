$(document).ready(function() {

    applyWhenElementExists(
        $('#initializing-block'),
        $('#initializing-modal'),
        [$('#init-step'), $('#init-message')],
        '/user/user-init-step',
        update_init,
        1000);

    $(".learn-more").click(function() {
       scrollToAnchor('learn-more');
    });
});


function applyWhenElementExists($block, $modal, $messages, url, myFunction, intervalTime) {
    // Control update of block / modal based on ajax calls

    var runningInterval = $block.data('interval');
    // if interval running, clear it
    if(runningInterval) {
        clearInterval(runningInterval);
        $block.removeData('interval');
    }
    // set interval
    var interval = setInterval(function() {
        if ($block.length > 0) {
            myFunction($block, $modal, $messages, url);
        }
    }, intervalTime);
    // store interval in object for latter clearing in myFunction
    $block.data('interval', interval);
}

function update_init($block,  $modal, $messages, url) {
    $.getJSON(url, function (json) {
        console.log(json);
        console.log($messages);
        if (json.done) {
            var libInterval = $block.data('interval');
            clearInterval(libInterval);
            $block.removeData('interval');
            $modal.hide()
                .siblings('#initializing-modal-backdrop')
                .hide();
            window.location.href = json.redirect;
        }
        else {
            $modal.show();
            $.each($messages, function (index, $sel) {
                console.log($sel);
                console.log(json.messages[index]);
                $sel.html(json.messages[index]);
            });
        }
    });
}


//function update_message(obj_up, obj_hide, url) {
//    $.getJSON(url, function (json) {
//        if (json.done) {
//            var libInterval = obj_up.data('interval');
//            clearInterval(libInterval);
//            obj_up.removeData('interval');
//            obj_hide.hide();
//            window.location.replace(json.url);
//        }
//        else {
//            obj_up.html(json.message);
//        }
//    });
//}

function scrollToAnchor(aid){
    var aTag = $("a[name='"+ aid +"']");
    $('html,body').animate({scrollTop: aTag.offset().top},'slow');
}