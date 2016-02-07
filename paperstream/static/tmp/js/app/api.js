define(['jquery'], function($) {

    var api = {};

    /**
     * Toggles the pinned state of the paper.
     *
     * @param id        (int) The paper id/pk
     * @param source    (string) The source page
     *
     * @returns Promise(resolve(pinned :bool), reject(error :string))
     */
    api.pin = function(id, source) {
        var d = $.Deferred();

        $.ajax({
            type: 'POST',
            url: '/user/paper/pin',
            data: {'pk': id, 'source': source}
        })
        .done(function(json) {
            if (json.hasOwnProperty('is_liked')) {
                d.resolve(true == json['is_liked']);
            } else {
                d.resolve(false);
            }
        })
        .fail(function() {
            d.reject('Pin request failed');
        });

        return d.promise();
    };

    /**
     * Toggles the banned state of the paper.
     *
     * @param id        (int) The paper id/pk
     * @param source    (string) The source page
     *
     * @returns Promise(resolve(banned :bool), reject(error :string))
     */
    api.ban = function(id, source) {
        var d = $.Deferred();

        $.ajax({
            type: 'POST',
            url: '/user/paper/ban',
            data: {'pk': id, 'source': source}
        })
        .done(function(json) {
            if (json.hasOwnProperty('is_ticked')) {
                d.resolve(true == json['is_ticked']);
            } else {
                d.resolve(false);
            }
        })
        .fail(function() {
            d.reject('Ban request failed');
        });

        return d.promise();
    };

    /**
     * Adds the paper to the user library.
     *
     * @param id    (int) The paper id/pk
     *
     * @returns Promise(resolve(added :bool), reject(error :string))
     */
    api.add = function(id) {
        var d = $.Deferred();

        $.ajax({
            type: 'POST',
            url: '/user/paper/add',
            data: {'pk': id}
        })
        .done(function(json) {
            if (json.hasOwnProperty('success')) {
                d.resolve(true == json['success']);
            } else {
                d.resolve(false);
            }
        })
        .fail(function() {
            d.reject('Add request failed');
        });

        return d.promise();
    };

    /**
     * Trashes the paper from the user library.
     *
     * @param id    (int) The paper id/pk
     *
     * @returns Promise(resolve(added :bool), reject(error :string))
     */
    api.trash = function(id) {
        var d = $.Deferred();

        $.ajax({
            type: 'POST',
            url: '/user/paper/trash',
            data: {'pk': id}
        })
        .done(function(json) {
            if (json.hasOwnProperty('success')) {
                d.resolve(true == json['success']);
            } else {
                d.resolve(false);
            }
        })
        .fail(function() {
            d.reject('Trash request failed');
        });

        return d.promise();
    };

    /**
     * Restores the paper into the user library.
     *
     * @param id    (int) The paper id/pk
     *
     * @returns Promise(resolve(added :bool), reject(error :string))
     */
    api.restore = function(id) {
        var d = $.Deferred();

        $.ajax({
            type: 'POST',
            url: '/user/paper/restore',
            data: {'pk': id}
        })
        .done(function(json) {
            if (json.hasOwnProperty('success')) {
                d.resolve(true == json['success']);
            } else {
                d.resolve(false);
            }
        })
        .fail(function() {
            d.reject('Trash request failed');
        });

        return d.promise();
    };

    return api;
});
