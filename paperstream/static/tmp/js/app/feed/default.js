define(['jquery', 'app/api', 'app/ui/list', 'bootstrap'], function($, Api) {

    $(function() {

        $('body')
            .on('etalia.publication.pin', function(e, data) {
                if (data.hasOwnProperty('is_liked')) {
                    $('.thumb[data-id=' + data.id + ']')
                        .find('.thumb-pin')
                        .toggleClass('active', data['is_liked']);
                }
            })
            .on('etalia.publication.ban', function(e, data) {
                if (data.hasOwnProperty('is_ticked') && data['is_ticked']) {
                    /*$('.thumb[data-id=' + data.id + ']')
                     .find('.thumb-ban')
                     .toggleClass('active', data['success']);*/

                    $('.thumb[data-id=' + data.id + ']').remove();
                }
            })
            .on('etalia.publication.trash', function(e, data) {
                if (data.hasOwnProperty('success') && data['success']) {
                    $('.thumb[data-id=' + data.id + ']').remove();
                }
            });
    });
});
