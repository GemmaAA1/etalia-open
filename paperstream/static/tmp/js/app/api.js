define(['jquery'], function($) {

    var $body, api = {
        debug: false
    };

    api.log = function() {
        if (this.debug) {
            console.log('[API] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }
    };

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
            api.log('pin success', data);
            data.id = id;
            $body.trigger('etalia.publication.pin', data);
        })
        .fail(function(xrh, status, error) {
            api.log('pin failure', xrh, status, error);
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
            api.log('ban success', data);
            data.id = id;
            $body.trigger('etalia.publication.ban', data);
        })
        .fail(function(xrh, status, error) {
            api.log('ban failure', xrh, status, error);
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
            api.log('add success', data);
            data.id = id;
            $body.trigger('etalia.publication.add', data);
        })
        .fail(function(xrh, status, error) {
            api.log('add failure', xrh, status, error);
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
            api.log('trash success', data);
            data.id = id;
            $body.trigger('etalia.publication.trash', data);
        })
        .fail(function(xrh, status, error) {
            api.log('trash failure', xrh, status, error);
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
            api.log('restore success', data);
            data.id = id;
            $body.trigger('etalia.publication.restore', data);
        })
        .fail(function(xrh, status, error) {
            api.log('restore failure', xrh, status, error);
        });
    };

    $(function() {
        $body = $('body');
    });

    return api;
});
