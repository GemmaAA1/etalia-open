define(['jquery'], function($) {

    var api = {}, $body;

    /**
     * Toggles the pinned state of the paper.
     *
     * @param id        (int) The paper id/pk
     * @param source    (string) The source page
     */
    api.pin = function(id, source) {
        source = source || '';
        $.ajax({
            type: 'POST',
            url: '/user/paper/pin',
            data: {'pk': id, 'source': source}
        })
        .done(function(data) {
            data.id = id;
            $body.trigger('etalia.publication.pin', data);
        })
        .fail(function() {
            console.log('Pin request failed');
        });
    };

    /**
     * Toggles the banned state of the paper.
     *
     * @param id        (int) The paper id/pk
     * @param source    (string) The source page
     */
    api.ban = function(id, source) {
        source = source || '';
        $.ajax({
            type: 'POST',
            url: '/user/paper/ban',
            data: {'pk': id, 'source': source}
        })
        .done(function(data) {
            data.id = id;
            $body.trigger('etalia.publication.ban', data);
        })
        .fail(function() {
            console.log('Ban request failed');
        });
    };

    /**
     * Adds the paper to the user library.
     *
     * @param id (int) The paper id/pk
     */
    api.add = function(id) {
        $.ajax({
            type: 'POST',
            url: '/user/paper/add',
            data: {'pk': id}
        })
        .done(function(data) {
            data.id = id;
            $body.trigger('etalia.publication.add', data);
        })
        .fail(function() {
            console.log('Add request failed');
        });
    };

    /**
     * Trashes the paper from the user library.
     *
     * @param id (int) The paper id/pk
     */
    api.trash = function(id) {
        $.ajax({
            type: 'POST',
            url: '/user/paper/trash',
            data: {'pk': id}
        })
        .done(function(data) {
            data.id = id;
            $body.trigger('etalia.publication.trash', data);
        })
        .fail(function() {
            console.log('Trash request failed');
        });
    };

    /**
     * Restores the paper into the user library.
     *
     * @param id (int) The paper id/pk
     */
    api.restore = function(id) {
        $.ajax({
            type: 'POST',
            url: '/user/paper/restore',
            data: {'pk': id}
        })
        .done(function(data) {
            data.id = id;
            $body.trigger('etalia.publication.restore', data);
        })
        .fail(function() {
            console.log('Trash request failed');
        });
    };

    $(function() {
        $body = $('body');
    });

    return api;
});
