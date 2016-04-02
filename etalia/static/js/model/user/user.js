define([
    'app/config',
    'backbone',
    'backbone-relational'
], function (Config, Backbone) {

    var UserModel = Backbone.RelationalModel.extend({
        urlRoot: Config.api_root + '/user/users',

        defaults: {
            email: null,
            first_name: null,
            last_name: null,
            photo_url: null
        },

        schema: {
            email: {type: 'Text', validators: ['required']},
            fisrt_name: {type: 'Text', validators: ['required']},
            last_name: {type: 'Text', validators: ['required']}
            //photo_url: {type: 'Text', validators: ['required']}
        }
    });

    /**
     * Returns the current (authenticated) user.
     *
     * @return UserModel
     *
     * @TODO fetch user id
     */
    UserModel.getCurrent = function() {
        var user = UserModel.find(21);
        if (!user) {
            user = new UserModel({id: 21});
            user.fetch();
        }
        return user;
    };

    return UserModel;
});
