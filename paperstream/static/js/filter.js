$(document).ready(function() {

    $('.like-img').on('click', function () {
        $(this).toggleClass('active');
        send_filter(undefined, undefined);
    });

    $('.checkbox').on('click', function () {
        send_filter();
    });

    //$('.journals .only').on('click', function (event) {
    //    $(this).closest('.side-panel-section').find('input').each(function () {
    //        $(this).prop('checked', false);
    //    });
    //    $(this).parents().siblings('input').prop('checked', true);
    //    send_filter('none', undefined);
    //    event.stopPropagation();
    //});
    //
    //$('.authors .only').on('click', function (event) {
    //    $(this).closest('.side-panel-section').find('input').each(function () {
    //        $(this).prop('checked', false);
    //    });
    //    $(this).parents().siblings('input').prop('checked', true);
    //    send_filter(undefined, 'none');
    //    event.stopPropagation();
    //});
    //
    //$('.journals .all').on('click', function () {
    //    $(this).closest('.side-panel-section').find('input').each(function () {
    //        $(this).prop('checked', true);
    //    });
    //    send_filter('all', undefined);
    //});
    //
    //$('.journals .none').on('click', function () {
    //    $(this).closest('.side-panel-section').find('input').each(function () {
    //        $(this).prop('checked', false);
    //    });
    //    send_filter('none', undefined);
    //});

    //$('.authors .all').on('click', function () {
    //    $(this).closest('.side-panel-section').find('input').each(function () {
    //        $(this).prop('checked', true);
    //    });
    //    send_filter(undefined, 'all');
    //});
    //
    //$('.authors .none').on('click', function () {
    //    $(this).closest('.side-panel-section').find('input').each(function () {
    //        $(this).prop('checked', false);
    //    });
    //    send_filter(undefined, 'none');
    //});

    //$('.journal, .author')
    //    .on('mouseenter', function () {
    //        $(this).find('.short').show();
    //        $(this).find('.long').hide();
    //    })
    //    .on('mouseleave', function () {
    //        $(this).find('.short').hide();
    //        $(this).find('.long').show();
    //    });

    // display first block at load
    $('ul').find('[data-block="1"]').show();

    //// show more block by block
    $('.show-more').on('click', function () {
         //get display block
        var last_vis_block = $(this).siblings('ul')
            .children('li:visible:last')
            .data('block');
        // hide show-more when reaching bottom
        var last_block = $(this).siblings('ul')
            .children('li:last')
            .data('block');
        var next = parseInt(last_vis_block) + 1;
        $(this).siblings('ul')
            .find("[data-block='" + next + "']")
            .show();

        if (last_block == next) {
            $(this).hide();
        }
         //show show-less
        $(this).siblings('.show-less').show();
    });

    // display show-less and bind onclick
    $('.show-less').on('click', function () {
        // get last vis block
        var last_vis_block = $(this).siblings('ul')
            .children('li:visible:last')
            .data('block');
        // hide block
        $(this).siblings('ul')
            .find("[data-block='" + last_vis_block + "']")
            .hide();
        // add show more if not visible
        $(this).siblings('.show-more').show();
        // remove show-less if reach top
        if (last_vis_block == "2") {
            $(this).hide();
        }
    });

});

//function send_filter(journals_flag, authors_flag, sorting_flag) {
function send_filter() {

    var json_data = {};
    json_data.source = 'filter';
    //json_data.journals_flag = journals_flag;
    //json_data.authors_flag = authors_flag;
    //json_data.sorting_flag = sorting_flag;
    json_data.like_flag = false;
    json_data.journals = [];
    json_data.authors = [];

    // build json object
    $('.journals').find('input').each(function () {
        if ($(this).is(':checked')) {
            //json_data['journals'].push([$(this).data('journal-title'), true]);
            json_data['journals'].push($(this).data('journal-title'));
        }
        //else {
        //    json_data['journals'].push([$(this).data('journal-title'), false]);
        //}
    });

    $('.authors').find('input').each(function () {
        if ($(this).is(':checked')) {
            //json_data['authors'].push([$(this).data('author-pk'), true]);
            json_data['authors'].push($(this).data('author-pk'));
        }
        //else {
        //    json_data['authors'].push([$(this).data('author-pk'), false]);
        //}
    });

    if ($('.like-img').hasClass('active')) {
        json_data.like_flag = true;
    }

    console.log(JSON.stringify(json_data));

    // ajax call and call back
    $.get(window.location.href, JSON.stringify(json_data), function (fragment) {
        $('.endless_data').html(fragment);
        bind_all();
        $.getScript( "https://d1bxh8uas1mnw7.cloudfront.net/assets/embed.js", function( data, textStatus, jqxhr ) {
                console.log('load performed');
        });
    })
}